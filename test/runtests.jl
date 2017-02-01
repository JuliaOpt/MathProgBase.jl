using GLPKMathProgInterface, Ipopt, ECOS

include("linprog.jl")
linprogtest(GLPKSolverLP())

include("mixintprog.jl")
mixintprogtest(GLPKSolverMIP())

include("quadprog.jl")
quadprogtest(IpoptSolver())
qpdualtest(IpoptSolver())
socptest(ECOSSolver())

include("nlp.jl")
nlptest(IpoptSolver())
nlptest_nohessian(IpoptSolver())
convexnlptest(IpoptSolver())
rosenbrocktest(IpoptSolver())

include("conicinterface.jl")
coniclineartest(ECOSSolver(), duals=true)
# Test conic fallback for LPs
coniclineartest(GLPKSolverLP(), duals=true)

include("linproginterface.jl")
linprogsolvertest(GLPKSolverLP())
linprogsolvertestextra(GLPKSolverLP())
# Test LP fallback for conics
linprogsolvertest(ECOSSolver())
