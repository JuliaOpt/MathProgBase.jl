# test default solvers

include("linprog.jl")
linprogtest()

include("mixintprog.jl")
mixintprogtest()

include("quadprog.jl")
quadprogtest()

# Test conic fallback for LPs
include("conicinterface.jl")
coniclineartest(MathProgBase.defaultLPsolver)
