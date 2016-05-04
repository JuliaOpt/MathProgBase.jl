#############################################################################
# MathProgBase.jl
#############################################################################
# test/conininterface.jl
# Test the MathProgBase.jl interface for a provided conic solver.
#############################################################################

using Base.Test
using MathProgBase

function coniclineartest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)

    # Problem LIN1 - all vars in nonneg cone
    # min -3x - 2y - 4z
    # st    x +  y +  z == 3 '4' -> z=2, x->2, obj -> -14
    #            y +  z == 2
    #       x>=0 y>=0 z>=0
    # Opt solution = -11
    # x = 1, y = 0, z = 2
    println("Problem LIN1")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
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
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] 3.0 tol
        @test_approx_eq_eps d[2] 1.0 tol
    end

    # Problem LIN1A - same as Problem LIN1, but with variable bounds
    #              as constraints instead
    # min   -3x - 2y - 4z
    # st  [3] - [x + y + z] ZERO
    #     [2] - [    y + z] ZERO
    #     [0] - [x        ] NONPOS   x>=0 -> 0<=x -> 0-x<=0
    #     [0] - [    y    ] NONPOS   y>=0 -> 0<=y -> 0-y<=0
    #     [0] - [       -z] NONNEG   z>=0 -> 0<=z -> 0-z<=0 -> 0+z>=0
    # Opt solution = -11
    # x = 1, y = 0, z = 2
    println("Problem LIN1A")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
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
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1]  3.0 tol
        @test_approx_eq_eps d[2]  1.0 tol
        @test_approx_eq_eps d[3]  0.0 tol
        @test_approx_eq_eps d[4] -2.0 tol
        @test_approx_eq_eps d[5]  0.0 tol
    end

    # Problem LIN2 - mixed free, nonneg, nonpos, zero, shuffled cones
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
    println("Problem LIN2")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
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
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] -7.0 tol
        @test_approx_eq_eps d[2] -2.0 tol
        @test_approx_eq_eps d[3]  4.0 tol
    end

    # Problem LIN2A - Problem LIN2 but with y,z variable bounds as constraints
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
    println("Problem LIN2A")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
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
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] -7.0 tol
        @test_approx_eq_eps d[2] -2.0 tol
        @test_approx_eq_eps d[3]  4.0 tol
        @test_approx_eq_eps d[4]  0.0 tol
        @test_approx_eq_eps d[5]  0.0 tol
    end

    # Problem LIN3 - Infeasible LP
    # min  0
    # s.t. x ≥ 1
    #      x ≤ -1
    # in conic form:
    # min 0
    # s.t. -1 + x ∈ R₊
    #       1 + x ∈ R₋
    println("Problem LIN3")

    b = [-1, 1]
    A = [-1, -1]''
    c = [0]
    constr_cones = [(:NonNeg,1:1),(:NonPos,2:2)]
    var_cones = [(:Free,1:1)]

    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Infeasible
    if duals
        y = MathProgBase.getdual(m)
        @test y[1] > 0
        @test y[2] < 0
        @test_approx_eq_eps A'y zeros(1) tol
        @test -dot(b,y) > 0
    end

    # Problem LIN4 - Infeasible LP
    # min  0
    # s.t. x ≥ 1
    #      x ≤ 0
    # in conic form:
    # min 0
    # s.t. -1 + x ∈ R₊
    #           x ∈ R₋
    println("Problem LIN4")

    b = [-1]
    A = [-1]''
    c = [0]
    constr_cones = [(:NonNeg,1:1)]
    var_cones = [(:NonPos,1:1)]

    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Infeasible
    if duals
        y = MathProgBase.getdual(m)
        @test y[1] > 0
        @test -dot(b,y) > 0
    end

end

function conicSOCtest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    # Problem SOC1
    # min 0x - 1y - 1z
    #  st  x            == 1
    #      x >= ||(y,z)||
    println("Problem SOC1")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
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
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] sqrt(2.0) tol
    end


    # Problem SOC1A - Problem SOC1 but in ECOS form
    # min      0x - 1y - 1z
    #  st [1]-[ x          ] ZERO
    #     [0]-[-x          ] SOC
    #     [0]-[     -y     ] SOC
    #     [0]-[          -z] SOC
    println("Problem SOC1A")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
    [ 0.0, -1.0, -1.0],
    [ 1.0   0.0   0.0;
     -1.0   0.0   0.0;
      0.0  -1.0   0.0;
      0.0   0.0  -1.0],
    [ 1.0, 0.0, 0.0, 0.0],
    Any[(:Zero,1),(:SOC,2:4)],
    [(:Free,1:3)])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) tol
    if duals
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] sqrt(2.0) tol
        @test_approx_eq_eps d[2] sqrt(2.0) tol
        @test_approx_eq_eps d[3] -1.0 tol
        @test_approx_eq_eps d[4] -1.0 tol
    end

    # Problem SOC2
    # min  x
    # s.t. y ≥ 1/√2
    #      x² + y² ≤ 1
    # in conic form:
    # min  x
    # s.t.  -1/√2 + y ∈ R₊
    #        1 - t ∈ {0}
    #      (t,x,y) ∈ SOC₃
    println("Problem SOC2")

    b = [-1/sqrt(2), 1, 0, 0, 0]
    A = [ 0 -1 0
          0 0 1
          0 0 -1
         -1 0 0
          0 -1 0 ]
    c = [ 1, 0, 0 ]
    constr_cones = [(:NonNeg,1:1),(:Zero,2:2),(:SOC,3:5)]
    var_cones = [(:Free,1:3)]

    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) -1/sqrt(2) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] -1/sqrt(2)  tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2]  1/sqrt(2)  tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3]  1.0  tol

    if duals
        y = MathProgBase.getdual(m)
        @test_approx_eq_eps -dot(b,y) -1/sqrt(2) tol
        @test_approx_eq_eps c + A'y zeros(3) tol
    end

    # Problem SOC2A
    # Same as above but with NonPos instead of NonNeg
    # min  x
    # s.t.  1/√2 - y ∈ R₋
    #        1 - t ∈ {0}
    #      (t,x,y) ∈ SOC₃
    println("Problem SOC2")

    b = [1/sqrt(2), 1, 0, 0, 0]
    A = [ 0 1 0
          0 0 1
          0 0 -1
          -1 0 0
          0 -1 0 ]
    c = [ 1, 0, 0 ]
    constr_cones = [(:NonPos,1:1),(:Zero,2:2),(:SOC,3:5)]
    var_cones = [(:Free,1:3)]

    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) -1/sqrt(2) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] -1/sqrt(2)  tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2]  1/sqrt(2)  tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3]  1.0  tol

    if duals
        y = MathProgBase.getdual(m)
        @test_approx_eq_eps -dot(b,y) -1/sqrt(2) tol
        @test_approx_eq_eps c + A'y zeros(3) tol
    end

    # Problem SOC3 - Infeasible
    # min 0
    # s.t. y ≥ 2
    #      x ≤ 1
    #      |y| ≤ x
    # in conic form:
    # min 0
    # s.t. -2 + y ∈ R₊
    #      -1 + x ∈ R₋
    #       (x,y) ∈ SOC₂
    println("Problem SOC3")
    b = [-2, -1]
    A = [0 -1; -1 0]
    c = [0,0]
    constr_cones = [(:NonNeg,1:1),(:NonPos,2:2)]
    var_cones = [(:SOC,1:2)]

    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Infeasible

    if duals
        y = MathProgBase.getdual(m)
        @test y[1] > 0
        @test y[2] < 0
        @test (A'y)[1] ≥ abs((A'y)[2])
        @test -dot(b,y) > 0
    end

end

function conicSOCINTtest(s::MathProgBase.AbstractMathProgSolver;tol=1e-6)
    # Problem SOCINT1
    # min 0x - 2y - 1z
    #  st  x            == 1
    #      x >= ||(y,z)||
    #      (y,z) binary
    println("Problem SOCINT1")
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m,
    [ 0.0, -2.0, -1.0],
    [ 1.0   0.0   0.0],
    [ 1.0],
    [(:Zero,1)],
    [(:SOC,1:3)])
    MathProgBase.setvartype!(m, [:Cont,:Bin,:Bin])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @show MathProgBase.getobjval(m)
    @show MathProgBase.getsolution(m)
    @test_approx_eq_eps MathProgBase.getobjval(m) -2 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3] 0.0 tol
end


function conicEXPtest(s::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    # Problem 4 - ExpPrimal
    # min x + y + z
    #  st  y e^(x/y) <= z, y > 0 (i.e (x, y, z) are in the exponential primal cone)
    #      x == 1
    #      y == 2

    println("Problem 4")
    m = MathProgBase.ConicModel(s)
    c = [1.0, 1.0, 1.0]
    A = [0.0 1.0 0.0;
         1.0 0.0 0.0]
    b = [2.0, 1.0]
    MathProgBase.loadproblem!(m, c, A, b,
    [(:Zero,1:2)],
    [(:ExpPrimal, 1:3)])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) (2*exp(1/2)+3) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 2.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3] 2*exp(1/2) tol

    if duals
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps -dot(b,d) MathProgBase.getobjval(m) tol
        u,v,w = c+A'd
        # should belong to the ExpDual cone
        @test u < 0
        @test -u*exp(v/w) ≤ exp(1)*w + tol
    end

    # Problem 4A - ExpPrimal
    # Same as previous, except we put :ExpPrimal on constr_cones
    println("Problem 4A")
    m = MathProgBase.ConicModel(s)
    c = [1.0, 1.0, 1.0]
    A = [0.0 1.0 0.0;
         1.0 0.0 0.0;
        -1.0 0.0 0.0;
         0.0 -1.0 0.0;
         0.0 0.0 -1.0]
    b = [2.0, 1.0, 0.0, 0.0, 0.0]
    MathProgBase.loadproblem!(m, c, A, b,
    [(:Zero,1:2), (:ExpPrimal, 3:5)],
    [(:Free, 1:3)])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test_approx_eq_eps MathProgBase.getobjval(m) (2*exp(1/2)+3) tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 2.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[3] 2*exp(1/2) tol

    if duals
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps -dot(b,d) MathProgBase.getobjval(m) tol
        @test_approx_eq_eps c+A'd zeros(3) tol
        u,v,w = d[3:5]
        # should belong to the ExpDual cone
        @test u < 0
        @test -u*exp(v/w) ≤ exp(1)*w + tol
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
    m = MathProgBase.ConicModel(s)
    c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
    A = [-eye(9); [0, 0, 0, 1, 0, 0, 0, 0, 0]']
    b = vcat(zeros(Int, 9), 1)
    MathProgBase.loadproblem!(m, c, A, b, [(:SDP, 1:9), (:Zero, 10:10)], [(:Free, 1:9)])
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
        d = MathProgBase.getdual(m)
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
    m = MathProgBase.ConicModel(s)
    c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
    A = [0, 0, 0, 1, 0, 0, 0, 0, 0]'
    b = [1]
    MathProgBase.loadproblem!(m, c, A, b, [(:Zero, 1:1)], [(:SDP, 1:9)])
    MathProgBase.optimize!(m)
    @test MathProgBase.status(m) == :Optimal
    @test is_symmetric(reshape(MathProgBase.getsolution(m), 3, 3))
    @test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0 tol
    @test_approx_eq_eps MathProgBase.getsolution(m)[4] 1.0 tol
    if duals
        d = MathProgBase.getdual(m)
        @test_approx_eq_eps d[1] 1 tol
    end

end
