using CTDE

typealias Time Float64

function within_house_model(params)
  cnt=params["count"]
  structure=ExplicitGSPN()

  add_place(structure, 's')

  birth=ConstExplicitTransition(
      (lm, user, when::Time)->begin
          (TransitionExponential(params["birth"](1-lm["s"]/params["carrying"]), when), Int[])
      end)
  add_transition(structure, 'b', birth,
      [TransitionRoute( 's', 's', )])

end
