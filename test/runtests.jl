using GLPKMathProgInterface, Ipopt, ECOS

include("linprog.jl")
linprogtest(GLPKSolverLP())

include("mixintprog.jl")
mixintprogtest(GLPKSolverMIP())

include("quadprog.jl")
quadprogtest(IpoptSolver())

include("conicinterface.jl")
coniclineartest(ECOSSolver(), duals=true)
# Test conic fallback for LPs
coniclineartest(GLPKSolverLP(), duals=true)

include("linproginterface.jl")
linprogsolvertest(GLPKSolverLP())
# Test LP fallback for conics
linprogsolvertest(ECOSSolver())
