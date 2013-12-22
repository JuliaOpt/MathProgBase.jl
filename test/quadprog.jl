using Base.Test
using MathProgBase

function quadprogtest(solver=MathProgBase.defaultQPsolver)
    println("Testing quadprog with solver ", string(typeof(solver)))

    sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],[4., 1.],Inf,-Inf,Inf,solver)
    @test sol.status == :Optimal
    @test_approx_eq sol.objval 130/70
    @test_approx_eq norm(sol.sol[1:3] - [0.5714285714285715,0.4285714285714285,0.8571428571428572]) 0.0

    println("Done")
end
