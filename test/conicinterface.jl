#############################################################################
# MathProgBase.jl
#############################################################################
# test/conininterface.jl
# Test the MathProgBase.jl interface for a provided conic solver.
#############################################################################

using Base.Test
using MathProgBase
using MathProgBase.MathProgSolverInterface

function coniclineartest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)

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
@test_approx_eq_eps MathProgBase.getobjval(m) -11 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -3.0 tol
  @test_approx_eq_eps d[2] -1.0 tol
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
@test_approx_eq_eps MathProgBase.getobjval(m) -11 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -3.0 tol
  @test_approx_eq_eps d[2] -1.0 tol
  @test_approx_eq_eps d[3]  0.0 tol
  @test_approx_eq_eps d[4] -2.0 tol
  @test_approx_eq_eps d[5]  0.0 tol
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
@test_approx_eq_eps MathProgBase.getobjval(m) -82 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] -4.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] -3.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 16.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[4]  0.0 tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1]  7.0 tol  # Check: change RHS to -3, x->-3, z->15
  @test_approx_eq_eps d[2]  2.0 tol  #        then obj +3+4=+7
  @test_approx_eq_eps d[3] -4.0 tol
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
@test_approx_eq_eps MathProgBase.getobjval(m) -82 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] -4.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] -3.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 16.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[4]  0.0 tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1]  7.0 tol
  @test_approx_eq_eps d[2]  2.0 tol
  @test_approx_eq_eps d[3] -4.0 tol
  @test_approx_eq_eps d[4]  0.0 tol
  @test_approx_eq_eps d[5]  0.0 tol
end

end

function conicSOCtest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
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
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -sqrt(2.0) tol
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
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] -sqrt(2.0) tol
  @test_approx_eq_eps d[2] -sqrt(2.0) tol
  @test_approx_eq_eps d[3] 1.0 tol
  @test_approx_eq_eps d[4] 1.0 tol
end

end

function conicEXPtest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
# Problem 4 - ExpPrimal
# min x + y + z
#  st  y e^(x/y) <= x, y > 0 (i.e (x, y, z) are in the exponential primal cone)
#      x == 1
#      y == 2

println("Problem 4")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
[1.0, 1.0, 1.0],
[0.0 1.0 0.0;
 1.0 0.0 0.0],
[2.0, 1.0],
[(:Zero,1:2)],
[(:ExpPrimal, 1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) 6.2974 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 2.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 3.29744 tol

if duals
  # cvx_begin
  #     variable x, y, z;
  #     dual variable a, b, c, d;
  #     minimize( x + y + z );
  #     subject to
  #         a : x == 1;
  #         b : y == 2
  #         { x, y, z } == exponential( 1 );
  # cvx_end
  # gives
  # a = 2.6488
  # b = 1.8243
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] 1.8243 tol
  @test_approx_eq_eps d[2] 2.6488 tol
end

# Problem 4A - ExpPrimal
# Same as previous, except we put :ExpPrimal on constr_cones
println("Problem 4A")
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
[1.0, 1.0, 1.0],
[0.0 1.0 0.0;
 1.0 0.0 0.0;
 -1.0 0.0 0.0;
 0.0 -1.0 0.0;
 0.0 0.0 -1.0],
[2.0, 1.0, 0.0, 0.0, 0.0],
[(:Zero,1:2), (:ExpPrimal, 3:5)],
[(:Free, 1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) 6.2974 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 2.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 3.29744 tol

if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] 1.8243 tol
  @test_approx_eq_eps d[2] 2.6488 tol
  # TODO: How do I calculate the dual of the exp constraint and check it with a non-SCS library?
  #
end
end


function conicSDPtest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)

function is_symmetric(A::Matrix)
  return all(A - A' .< 1e-4)
end

# Problem 5 - SDP
# min y[1, 2]
#  st y[2, 1] == 1
#     y in SDP cone
# If symmetricity constraint is working, y[1, 2] will be 1 else unbounded
m = MathProgBase.model(s)
c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
A = -eye(9)
A = [A; [0, 0, 0, 1, 0, 0, 0, 0, 0]']
b = zeros(size(A, 1), 1)
b[10] = 1
MathProgBase.loadconicproblem!(m, c, A, b, [(:SDP, 1:9), (:Zero, 10:10)], [(:Free, 1:9)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test is_symmetric(reshape(MathProgBase.getsolution(m), 3, 3))
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[4] 1.0 tol
if duals
  # cvx_begin
  #  variable y(3,3)
  #  dual variables a b;
  #  minimize y(1,2)
  #  subject to
  #  a : y(2,1) == 1;
  #  b : y == semidefinite(3)
  # cvx_end
  # gives
  # a = 1
  # b = [0 1 0; -1 0 0; 0 0 0]
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] 0 tol
  @test_approx_eq_eps d[2] -1 tol;
  @test_approx_eq_eps d[3] 0 tol;
  @test_approx_eq_eps d[4] 1 tol;
  @test_approx_eq_eps d[5] 0 tol;
  @test_approx_eq_eps d[6] 0 tol;
  @test_approx_eq_eps d[7] 0 tol;
  @test_approx_eq_eps d[8] 0 tol;
  @test_approx_eq_eps d[9] 0 tol;
  @test_approx_eq_eps d[10] 1 tol
end

# Problem 5A - SDP
# Same as problem 5, except we enforce :SDP on the var_cone
m = MathProgBase.model(s)
c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
A = [0, 0, 0, 1, 0, 0, 0, 0, 0]'
b = [1]
MathProgBase.loadconicproblem!(m, c, A, b, [(:Zero, 1:1)], [(:SDP, 1:9)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test is_symmetric(reshape(MathProgBase.getsolution(m), 3, 3))
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0 tol
@test_approx_eq_eps MathProgBase.getsolution(m)[4] 1.0 tol
if duals
  d = MathProgBase.getconicdual(m)
  @test_approx_eq_eps d[1] 1 tol
end

end
