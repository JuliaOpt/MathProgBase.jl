using Base.Test
using MathProgBase

solver = MIPSolver(:Cbc,LogLevel=1)

# integer knapsack problem
sol = mixintprog(-[5,3,2,7,4],Float64[2 8 4 2 5],'<',10,'I',0,1,solver)
@test sol.status == :Optimal
@test_approx_eq sol.objval -16.0
@test_approx_eq norm(sol.sol[1:5] - [1.0, 0.0, 0.0, 1.0, 1.0]) 0.0
