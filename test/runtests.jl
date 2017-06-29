using GLPKMathProgInterface, Ipopt, ECOS

const MPB = MathProgBase

lp_solver = GLPKSolverLP()
ip_solver = GLPKSolverMIP()
conic_solver = ECOSSolver(verbose=false)
nlp_solver = IpoptSolver(print_level=0)

include("linprog.jl")
linprogtest(lp_solver)

include("mixintprog.jl")
mixintprogtest(ip_solver)

include("quadprog.jl")
quadprogtest(nlp_solver)
qpdualtest(nlp_solver)
socptest(conic_solver)

include("nlp.jl")
nlptest(nlp_solver)
nlptest_nohessian(nlp_solver)
convexnlptest(nlp_solver)
rosenbrocktest(nlp_solver)

include("conicinterface.jl")
coniclineartest(conic_solver, duals=true)
# Test conic fallback for LPs
coniclineartest(lp_solver, duals=true)

include("linproginterface.jl")
linprogsolvertest(lp_solver)
linprogsolvertestextra(lp_solver)
# Test LP fallback for conics
linprogsolvertest(conic_solver)
