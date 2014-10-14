# test default solvers

include("linprog.jl")
linprogtest()

include("mixintprog.jl")
mixintprogtest()

include("quadprog.jl")
quadprogtest()


# Test fallback for loadconicproblem! on an LP
# Running tests for MathProgBase requires Clp anyway, so we'll
# use that as our solver
using Base.Test
using MathProgBase

function test_conic_fallback()
    m = MathProgBase.model(MathProgBase.defaultLPsolver)
    MathProgBase.loadconicproblem!(m, 
                        [-3.0, -2.0, -4.0],
                        [ 1.0   1.0   1.0;
                          0.0   1.0   1.0],
                        [ 3.0,  2.0],
                        [(:Zero,1:2)],
                        [(:NonNeg, 1:3)])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) -11 1e-6
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-6
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 1e-6
    @test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 1e-6
end
test_conic_fallback()