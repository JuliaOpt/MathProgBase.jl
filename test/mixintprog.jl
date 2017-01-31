using Base.Test
using MathProgBase

function mixintprogtest(solver)
    println("Testing mixintprog with solver ", string(typeof(solver)))

    # integer knapsack problem
    sol = mixintprog(-[5,3,2,7,4],Float64[2 8 4 2 5],'<',10,:Int,0,1,solver)
    @test sol.status == :Optimal
    @test isapprox(sol.objval, -16.0, atol=1e-6)
    @test isapprox(norm(sol.sol - [1.0, 0.0, 0.0, 1.0, 1.0]), 0.0, atol=1e-4)

    println("Done")
end
