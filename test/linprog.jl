using Base.Test
using MathProgBase

function linprogtest(solver=MathProgBase.defaultLPsolver; objtol = 1e-7, primaltol = 1e-6)
    println("Testing linprog with solver ", string(typeof(solver)))
    # min -x
    # s.t. 2x + y <= 1.5
    # x,y >= 0
    # solution is (0.75,0) with objval -0.75

    sol = linprog([-1,0],[2 1],'<',1.5,solver)
    @test sol.status == :Optimal
    @test_approx_eq_eps sol.objval -0.75 objtol
    @test_approx_eq_eps norm(sol.sol - [0.75,0.0]) 0 primaltol

    sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
    @test sol.status == :Optimal
    @test_approx_eq_eps sol.objval -0.75 objtol
    @test_approx_eq_eps norm(sol.sol - [0.75,0.0]) 0 primaltol

    # test infeasible problem:
    # min x
    # s.t. 2x+y <= -1
    # x,y >= 0
    sol = linprog([1,0],[2 1],'<',-1,solver)
    @test sol.status == :Infeasible

    r = sol.attrs[:infeasibilityray][1] 
    @test_approx_eq r/abs(r) -1.0

    # test unbounded problem:
    # min -x-y
    # s.t. -x+2y <= 0
    # x,y >= 0
    sol = linprog([-1,-1],[-1 2],'<',[0],solver)
    @test sol.status == :Unbounded

    # unbounded problem with unique ray:
    # min -x-y
    # s.t. x-y == 0
    # x,y >= 0
    sol = linprog([-1,-1],[1 -1],'=',0,solver)
    @test sol.status == :Unbounded
    @test sol.attrs[:unboundedray][1] > 1e-7
    @test_approx_eq sol.attrs[:unboundedray][1] sol.attrs[:unboundedray][2]


    println("Passed")
end



