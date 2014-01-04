# test default solvers

include("linprog.jl")
linprogtest()

include("mixintprog.jl")
mixintprogtest()

include("quadprog.jl")
quadprogtest()
