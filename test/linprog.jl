using Base.Test
using MathProgBase

function linprogtest(solver=MathProgBase.defaultLPSolver)
    println("Testing linprog with solver ", string(typeof(solver)))
    # min -x
    # s.t. 2x + y <= 1.5
    # x,y >= 0
    # solution is (0.75,0) with objval -0.75

    sol = linprog([-1,0],[2 1],'<',1.5,solver)
    @test sol.status == :Optimal
    @test_approx_eq sol.objval -0.75
    @test_approx_eq norm(sol.sol - [0.75,0.0]) 0 

    sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
    @test sol.status == :Optimal
    @test_approx_eq sol.objval -0.75
    @test_approx_eq norm(sol.sol - [0.75,0.0]) 0 

    # test infeasible problem:
    # min x
    # s.t. 2x+y <= -1
    # x,y >= 0
    sol = linprog([1,0],[2 1],'<',-1,solver)
    @test sol.status == :Infeasible

    # test unbounded problem:
    # min -x-y
    # s.t. -x+2y <= 0
    # x,y >= 0
    sol = linprog([-1,-1],[-1 2],'<',[0],solver)
    @test sol.status == :Unbounded

    println("Passed")
end



