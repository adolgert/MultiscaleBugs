using CTDE

typealias Time Float64

function hop_skip_graph(params, propensity_matrix)
    cnt=size(propensity_matrix)[1]
    structure=ExplicitGSPN()

    for ind_place_idx in 1:cnt
        add_place(structure, (ind_place_idx, 's'))
        add_place(structure, (ind_place_idx, 'i'))
        add_place(structure, (ind_place_idx, 'r'))
    end

    for alone_idx in 1:cnt
        recover=ConstExplicitTransition(
            (lm, user, when::Time)->begin
                (TransitionExponential(params["gamma"], when), Int[])
            end)
        add_transition(structure, (alone_idx, alone_idx, 'r'), recover,
            [TransitionRoute( (alone_idx, 'i'), (alone_idx, 'r'), 1),],
            [])
    end

    for spread_idx in 1:cnt
        for neighbor_idx in 1:cnt
            if neighbor_idx!=spread_idx
                beta=propensity_matrix[spread_idx, neighbor_idx]
                infect=ConstExplicitTransition(
                    (lm, user, when::Time)->begin
                        (TransitionExponential(beta, when), Int[])
                end)
                (i, j)=(spread_idx, neighbor_idx)
                add_transition(structure, (i, j, 'b'), infect,
                    [TransitionRoute((i, 'i'), (i, 'i'), 1),
                     TransitionRoute((j, 's'), (j, 'i'), 1)],
                     [])
            end
        end
    end
    structure
end


function initialize_marking(model, cnt, first)
    for individual_idx in 1:cnt
        if individual_idx!=first
            add_tokens(model, (individual_idx, 's'), 1)
        else
            add_tokens(model, (individual_idx, 'i'), 1)
        end
    end
end


type DiseaseObserver
    sir::Array{Time, 2}
    infected::Int
    DiseaseObserver(cnt)=new(zeros(Time, 2, cnt), 1)
end


function observe(eo::DiseaseObserver, state)
    last_fired=state.last_fired
    if last_fired[3]=='r'
        eo.sir[2, last_fired[1]]=state.current_time
        eo.infected-=1
    elseif last_fired[3]=='b'
        eo.sir[1, last_fired[2]]=state.current_time
        eo.infected+=1
    else
        error("Unknown transition")
    end
    eo.infected>0
end


function single_hop_skip(params, propensity_matrix, myseed)
    rng=MersenneTwister(myseed)
    cnt=size(propensity_matrix)[1]
    state=TokenState(int_marking(), params)
    structure=hop_skip_graph(params, propensity_matrix)
    model=GSPNModel(structure, state)

    sampling=NextReactionHazards()
    observer=DiseaseObserver(cnt)
    first=rand(1:cnt)
    initialize_marking(model, cnt, first)
    run_steps(model, sampling, s->observe(observer, s), rng)
    observer.sir
end

