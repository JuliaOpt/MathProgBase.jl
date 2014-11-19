#############################################################################
# MathProgBase.jl
#############################################################################
# test/conininterface.jl
# Test the MathProgBase.jl interface for a provided conic solver.
#############################################################################

using Base.Test
using MathProgBase
using MathProgBase.MathProgSolverInterface

function coniclineartest(s::MathProgBase.AbstractMathProgSolver;duals=false)

# Problem 1 - all vars in nonneg cone
# min -3x - 2y - 4z
# st    x +  y +  z == 3 '4' -> z=2, x->2, obj -> -14
#            y +  z == 2
#       x>=0 y>=0 z>=0
# Opt solution = -11
# x = 1, y = 0, z = 2
println("Problem 1")
m = MathProgBase.model(s)
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
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -3.0 1e-6
  @test_approx_eq_eps d[2] -1.0 1e-6
end

# Problem 1A - same as Problem 1, but with variable bounds
#              as constraints instead
# min   -3x - 2y - 4z
# st  [3] - [x + y + z] ZERO
#     [2] - [    y + z] ZERO
#     [0] - [x        ] NONPOS   x>=0 -> 0<=x -> 0-x<=0
#     [0] - [    y    ] NONPOS   y>=0 -> 0<=y -> 0-y<=0
#     [0] - [       -z] NONNEG   z>=0 -> 0<=z -> 0-z<=0 -> 0+z>=0
# Opt solution = -11
# x = 1, y = 0, z = 2
println("Problem 1A")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m, 
                    [-3.0, -2.0, -4.0],
                    [ 1.0   1.0   1.0;
                      0.0   1.0   1.0;
                      1.0   0.0   0.0;
                      0.0   1.0   0.0;
                      0.0   0.0  -1.0],
                    [ 3.0,  2.0,  0.0,  0.0,  0.0],
                    [(:Zero,1:2),(:NonPos,3:4),(:NonNeg,5)],
                    [(:Free, 1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -11 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 1e-6
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -3.0 1e-6
  @test_approx_eq_eps d[2] -1.0 1e-6
  @test_approx_eq_eps d[3]  0.0 1e-6
  @test_approx_eq_eps d[4] -2.0 1e-6
  @test_approx_eq_eps d[5]  0.0 1e-6
end

# Problem 2 - mixed free, nonneg, nonpos, zero, shuffled cones
# min  3x + 2y - 4z + 0s
# st    x           -  s  == -4    (i.e. x >= -4)
#            y            == -3
#       x      +  z       == 12
#       x free
#       y <= 0
#       z >= 0
#       s zero
# Opt solution = -82
# x = -4, y = -3, z = 16, s == 0
println("Problem 2")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 3.0,  2.0, -4.0,  0.0],
                    [ 1.0   0.0   0.0  -1.0;
                      0.0   1.0   0.0   0.0;
                      1.0   0.0   1.0   0.0],
                    [-4.0, -3.0, 12.0],
                    [(:Zero,1:3)],
                    [(:Free,1), (:NonNeg,3), (:Zero,4), (:NonPos,2)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -82 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] -4.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] -3.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 16.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[4]  0.0 1e-6
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1]  7.0 1e-6  # Check: change RHS to -3, x->-3, z->15
  @test_approx_eq_eps d[2]  2.0 1e-6  #        then obj +3+4=+7
  @test_approx_eq_eps d[3] -4.0 1e-6
end

# Problem 2A - Problem 2 but with y,z variable bounds as constraints
# min       3x + 2y - 4z + 0s
# st  [-4]-[ x           -  s] ZERO    (i.e. x >= -4)
#     [-3]-[      y          ] ZERO
#     [12]-[ x      +  z     ] ZERO
#     [ 0]-[      y          ] NONEG
#     [ 0]-[           z     ] NONPOS
#       x free
#       s zero
# Opt solution = -82
# x = -4, y = -3, z = 16, s == 0
println("Problem 2A")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 3.0,  2.0, -4.0,  0.0],
                    [ 1.0   0.0   0.0  -1.0;
                      0.0   1.0   0.0   0.0;
                      1.0   0.0   1.0   0.0;
                      0.0   1.0   0.0   0.0;
                      0.0   0.0   1.0   0.0],
                    [-4.0, -3.0, 12.0,  0.0,  0.0],
                    [(:Zero,1:3),(:NonNeg,4),(:NonPos,5)],
                    [(:Free,1), (:NonNeg,3), (:Zero,4), (:NonPos,2)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -82 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] -4.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] -3.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 16.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[4]  0.0 1e-6
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1]  7.0 1e-6
  @test_approx_eq_eps d[2]  2.0 1e-6
  @test_approx_eq_eps d[3] -4.0 1e-6
  @test_approx_eq_eps d[4]  0.0 1e-6
  @test_approx_eq_eps d[5]  0.0 1e-6
end

end

function conicSOCtest(s::MathProgBase.AbstractMathProgSolver;duals=false)
# Problem 3 - SOC
# min 0x - 1y - 1z
#  st  x            == 1
#      x >= ||(y,z)||
println("Problem 3")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 0.0, -1.0, -1.0],
                    [ 1.0   0.0   0.0],
                    [ 1.0],
                    [(:Zero,1)],
                    [(:SOC,1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) 1e-6
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -sqrt(2.0) 1e-6
end


# Problem 3A - Problem 3 but in ECOS form
# min      0x - 1y - 1z
#  st [1]-[ x          ] ZERO 
#     [0]-[-x          ] SOC
#     [0]-[     -y     ] SOC
#     [0]-[          -z] SOC
println("Problem 3A")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 0.0, -1.0, -1.0],
                    [ 1.0   0.0   0.0;
                     -1.0   0.0   0.0;
                      0.0  -1.0   0.0;
                      0.0   0.0  -1.0],
                    [ 1.0, 0.0, 0.0, 0.0],
                    {(:Zero,1),(:SOC,2:4)},
                    [(:Free,1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) 1e-6
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) 1e-6
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -sqrt(2.0) 1e-6
  @test_approx_eq_eps d[2] -sqrt(2.0) 1e-6
  @test_approx_eq_eps d[3] 1.0 1e-6
  @test_approx_eq_eps d[4] 1.0 1e-6
end

end
