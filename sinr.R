# Simulate susceptible-infected-notified-recovered with Gillespie.

# This is the library from which we stole the Gillespie Direct method.
#library(smfsb)


# Complete spatial randomness point process
csr_process <- function(N) {
  #locations uniform on [0,1]x[0,1] grid
  location<-matrix(c(NA,NA),nrow=N,ncol=2) 
  for (i in 1:N) {location[i,]<-c(runif(1),runif(1))}
  location
}


# Given a list of x-y pairs, return a distance matrix.
distance_matrix <- function(location) {
  N=dim(location)[[1]]
  distance<-matrix(NA, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      distance[i, j]=sqrt((location[i,1]-location[j,1])^2+
        (location[i,2]-location[j,2])^2)
    }
  }
  distance
}


# Build a continuous-time SIR model given a distance matrix.
# This is a chemical-kinetics style of specification.
sir_model <- function(distance, parameters) {
  cat("sir_model\n")
  N<-dim(distance)[[1]]
  reaction_cnt<-N*N + N
  cat("N", N, "\n")

  # The species are stored in a vector 3N long
  # for N susceptible, N infected, N recovered.
  lhs=matrix(data=rep(0.0, reaction_cnt*3*N), nrow=reaction_cnt, ncol=3*N)
  rhs=matrix(data=rep(0.0, reaction_cnt*3*N), nrow=reaction_cnt, ncol=3*N)
  rates=vector("numeric", reaction_cnt)

  cat("Starting infections\n")
  infection_hazard=parameters$beta*ifelse(distance<0.3, 1, 0.1)
  reaction_idx=1
  for (inf_idx in 1:N) {
    for (sus_idx in 1:N) {

      # S + I -> I + I
      lhs[reaction_idx, inf_idx+N]=1
      lhs[reaction_idx, sus_idx]=1
      rhs[reaction_idx, inf_idx+N]=1
      rhs[reaction_idx, sus_idx+N]=1
      rates[reaction_idx]=infection_hazard[inf_idx, sus_idx]
      reaction_idx<-reaction_idx+1
    }
  }

  cat("recovery", reaction_idx, "\n")

  for (recover_idx in 1:N) {

    # I -> R
    lhs[reaction_idx, recover_idx+N]=1
    rhs[reaction_idx, recover_idx+2*N]=1
    rates[reaction_idx]<-parameters$gamma
    reaction_idx<-reaction_idx+1
  }

  print("created stoichiometry")
  # This creates a function which remembers local variables.
  # It returns a current rate for every infection and recovery.
  rate_function <- function(state, time) {
    current_rate=rep(0.0, N*N + N)
    reaction_idx=1
    for (inf_idx in 1:N) {
      if (state[inf_idx+N]>0) {
        for (sus_idx in 1:N) {
          if (state[sus_idx]>0) {
            current_rate[reaction_idx]=state[sus_idx]*rates[reaction_idx]
          }
          reaction_idx<-reaction_idx+1
        }
      } else {
        reaction_idx<-reaction_idx+N
      }
    }

    for (recover_idx in 1:N) {
      if (state[recover_idx+N]>0) {
        current_rate[reaction_idx]=rates[reaction_idx]
      }
      reaction_idx<-reaction_idx+1
    }
    current_rate
  }

  list(LHS=lhs, RHS=rhs, RATE=rate_function)
}



spatial_sir<-function(parameters, location_cnt) {
  distance<-distance_matrix(csr_process(location_cnt))
  sir_model(distance, parameters)
}



# Copied from gillespie() in smfsb package.
# This is a straight-up Gillespie Direct method.
# Modified to handle the possibility there are no more
# reactions at some step.
gillespie_copy <- function(N, n, ...)
{
  tt = 0
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  tvec = vector("numeric",n)
  xmat = matrix(nrow=n+1, ncol=u)
  cat("x", length(x), "\n")
  xmat[1,] = x
  for (i in 1:n) {
    h = N$h(x, tt, ...)
    total_propensity=sum(h)
    if (total_propensity>0){
        tt = tt+rexp(1,total_propensity)
        j = sample(v,1,prob=h)
        x = x+S[,j]
    }
    tvec[i] = tt
    xmat[i+1,] = x
  }
  return(list(t=tvec, x=xmat))
}


impute_observed_population<-function(trajectory, parameters) {
  t=trajectory$t
  x=trajectory$x

  time_cnt=length(t)
  compartment_cnt=dim(x)[[2]]
  individual_cnt=compartment_cnt/3

  observed=rep(0, individual_cnt)

  # Given a time since infection, this estimates how many bugs are observed.
  # The observation is stochastic, but the growth of bugs is not.
  observed_logistic<-function(time_delta) {
    ert=exp(parameters$growth * time_delta)
    population<-parameters$carrying * ert/(parameters$carrying+(ert-1))
    intpop=round(population)
    # observation_rate is the binomial chance of seeing a bug.
    observed<-rbinom(n=1, size=intpop, prob=parameters$observation_rate)
  }

  for (ind_idx in 1:individual_cnt) {
    infection_time=-1.0
    observation_time=-1.0
    for (time_idx in 1:time_cnt) {
      if (x[time_idx, ind_idx+individual_cnt]==1) {
        infection_time=t[time_idx]
        break
      }
    }
    for (time_idx in 1:time_cnt) {
      if (x[time_idx, ind_idx+2*individual_cnt]==1) {
        observation_time=t[time_idx]
        break
      }
    }
    if (infection_time>0 && observation_time>0) {
      observed[ind_idx]=observed_logistic(observation_time-infection_time)
    } else {
      observed[ind_idx]=0
    }
  }
  observed
}



run_one<-function(parameters, N) {
  model<-spatial_sir(parameters, N)
  marking<-vector("numeric", 3*N)
  marking[1:N]<-rep(1,N)
  marking[(N+1):(3*N)]<-rep(0, 2*N)
  initial_infective<-sample(1:N, 1)
  marking[initial_infective]=0
  marking[initial_infective+N]=1

  system<-list(M=marking, Pre=model$LHS, Post=model$RHS, h=model$RATE)
  cat("marking ", length(marking), "\n")
  print(dim(marking))
  cat("LHS", dim(model$LHS), "\n")
  cat("LHS row", nrow(model$LHS), "LHS ncol", ncol(model$LHS), "\n")
  print(dim(model$LHS))
  trajectory<-gillespie_copy(system, 2*N)

  observed<-impute_observed_population(trajectory, parameters)
  observed
}


set.seed(33333, kind=NULL, normal.kind=NULL)
run_one(list(beta=0.03, gamma=0.02, carrying=1000, growth=1/.9,
    observation_rate=0.5), 100)

