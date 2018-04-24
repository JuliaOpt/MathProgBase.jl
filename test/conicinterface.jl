#############################################################################
# MathProgBase.jl
#############################################################################
# test/conininterface.jl
# Test the MathProgBase.jl interface for a provided conic solver.
#############################################################################

using Compat.Test
using Compat.LinearAlgebra
using Compat.SparseArrays
using MathProgBase

function coniclineartest(solver::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    @testset "Testing linear problems through conic interface with $solver" begin
        # Problem LIN1 - all vars in nonneg cone
        # min -3x - 2y - 4z
        # st    x +  y +  z == 3 '4' -> z=2, x->2, obj -> -14
        #            y +  z == 2
        #       x>=0 y>=0 z>=0
        # Opt solution = -11
        # x = 1, y = 0, z = 2
        @testset "LIN1" begin
            c = [-3.0, -2.0, -4.0]
            A = [ 1.0   1.0   1.0;
                  0.0   1.0   1.0]
            b = [ 3.0,  2.0]
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, [(:Zero,1:2)], [(:NonNeg, 1:3)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -11, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 0.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 2.0, atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], 3.0, atol=tol)
                @test isapprox(d[2], 1.0, atol=tol)
                vardual = c + A'd
                var = MathProgBase.getvardual(m)
                @test isapprox(vardual, var, atol=tol)
            end
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
        @testset "LIN1A" begin
            c = [-3.0, -2.0, -4.0]
            A = [ 1.0   1.0   1.0;
                  0.0   1.0   1.0;
                  1.0   0.0   0.0;
                  0.0   1.0   0.0;
                  0.0   0.0  -1.0]
            b = [ 3.0,  2.0,  0.0,  0.0,  0.0]
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,c, A, b, [(:Zero,1:2),(:NonPos,3:4),(:NonNeg,5)],
            [(:Free, 1:3)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -11, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 0.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 2.0, atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], 3.0, atol=tol)
                @test isapprox(d[2], 1.0, atol=tol)
                @test isapprox(d[3], 0.0, atol=tol)
                @test isapprox(d[4], -2.0, atol=tol)
                @test isapprox(d[5], 0.0, atol=tol)
                vardual = c + A'd
                var = MathProgBase.getvardual(m)
                @test isapprox(vardual, var, atol=tol)
            end
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
        @testset "LIN2" begin
            m = MathProgBase.ConicModel(solver)
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
            @test isapprox(MathProgBase.getobjval(m), -82, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], -4.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], -3.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 16.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[4], 0.0, atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], -7.0, atol=tol)
                @test isapprox(d[2], -2.0, atol=tol)
                @test isapprox(d[3], 4.0, atol=tol)
            end
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
        @testset "LIN2A" begin
            m = MathProgBase.ConicModel(solver)
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
            @test isapprox(MathProgBase.getobjval(m), -82, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], -4.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], -3.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 16.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[4], 0.0, atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], -7.0, atol=tol)
                @test isapprox(d[2], -2.0, atol=tol)
                @test isapprox(d[3], 4.0, atol=tol)
                @test isapprox(d[4], 0.0, atol=tol)
                @test isapprox(d[5], 0.0, atol=tol)
            end
        end

        # Problem LIN3 - Infeasible LP
        # min  0
        # s.t. x ≥ 1
        #      x ≤ -1
        # in conic form:
        # min 0
        # s.t. -1 + x ∈ R₊
        #       1 + x ∈ R₋
        @testset "LIN3 infeasible" begin
            b = [-1, 1]
            A = reshape([-1, -1], 2, 1)
            c = [0]
            constr_cones = [(:NonNeg,1:1),(:NonPos,2:2)]
            var_cones = [(:Free,1:1)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Infeasible
            if duals
                y = MathProgBase.getdual(m)
                @test y[1] > 0
                @test y[2] < 0
                @test isapprox((A'y)[1], 0.0, atol=tol)
                @test -dot(b,y) > 0
            end
        end

        # Problem LIN4 - Infeasible LP
        # min  0
        # s.t. x ≥ 1
        #      x ≤ 0
        # in conic form:
        # min 0
        # s.t. -1 + x ∈ R₊
        #           x ∈ R₋
        @testset "LIN4 infeasible" begin
            b = [-1]
            A = reshape([-1], 1, 1)
            c = [0]
            constr_cones = [(:NonNeg,1:1)]
            var_cones = [(:NonPos,1:1)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Infeasible
            if duals
                y = MathProgBase.getdual(m)
                @test y[1] > 0
                @test -dot(b,y) > 0
            end
        end
    end
end

function conicSOCtest(solver::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    @testset "Testing SOC problems through conic interface with $solver" begin
        # Problem SOC1
        # min 0x - 1y - 1z
        #  st  x            == 1
        #      x >= ||(y,z)||
        @testset "SOC1" begin
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,
            [ 0.0, -1.0, -1.0],
            [ 1.0   0.0   0.0],
            [ 1.0],
            [(:Zero,1)],
            [(:SOC,1:3)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1.0/sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0/sqrt(2.0), atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], sqrt(2.0), atol=tol)
            end
        end


        # Problem SOC1A - Problem SOC1 but in ECOS form
        # min      0x - 1y - 1z
        #  st [1]-[ x          ] ZERO
        #     [0]-[-x          ] SOC
        #     [0]-[     -y     ] SOC
        #     [0]-[          -z] SOC
        @testset "SOC1A" begin
            m = MathProgBase.ConicModel(solver)
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
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1.0/sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0/sqrt(2.0), atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], sqrt(2.0), atol=tol)
                @test isapprox(d[2], sqrt(2.0), atol=tol)
                @test isapprox(d[3], -1.0, atol=tol)
                @test isapprox(d[4], -1.0, atol=tol)
            end
            # Permute indices - used to cause Mosek failure
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,
            [ 0.0, -1.0, -1.0],
            [ 1.0   0.0   0.0;
              0.0   0.0  -1.0;
              0.0  -1.0   0.0;
             -1.0   0.0   0.0],
            [ 1.0, 0.0, 0.0, 0.0],
            Any[(:Zero,1),(:SOC,[4,3,2])],
            [(:Free,1:3)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1.0/sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0/sqrt(2.0), atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(d[1], sqrt(2.0), atol=tol)
                @test isapprox(d[4], sqrt(2.0), atol=tol)
                @test isapprox(d[3], -1.0, atol=tol)
                @test isapprox(d[2], -1.0, atol=tol)
            end

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
        @testset "SOC2" begin
            b = [-1/sqrt(2), 1, 0, 0, 0]
            A = [ 0 -1 0
                  0 0 1
                  0 0 -1
                 -1 0 0
                  0 -1 0 ]
            c = [ 1, 0, 0 ]
            constr_cones = [(:NonNeg,1:1),(:Zero,2:2),(:SOC,3:5)]
            var_cones = [(:Free,1:3)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], -1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0, atol=tol)

            if duals
                y = MathProgBase.getdual(m)
                @test isapprox(-dot(b,y), -1/sqrt(2), atol=tol)
                @test isapprox(norm(c + A'y), 0.0, atol=tol)
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s), 0.0, atol=tol)
            end
        end

        # Problem SOC2A
        # Same as above but with NonPos instead of NonNeg
        # min  x
        # s.t.  1/√2 - y ∈ R₋
        #        1 - t ∈ {0}
        #      (t,x,y) ∈ SOC₃
        @testset "SOC2A" begin
            b = [1/sqrt(2), 1, 0, 0, 0]
            A = [ 0 1 0
                  0 0 1
                  0 0 -1
                  -1 0 0
                  0 -1 0 ]
            c = [ 1, 0, 0 ]
            constr_cones = [(:NonPos,1:1),(:Zero,2:2),(:SOC,3:5)]
            var_cones = [(:Free,1:3)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], -1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1/sqrt(2), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0, atol=tol)

            if duals
                y = MathProgBase.getdual(m)
                @test isapprox(-dot(b,y), -1/sqrt(2), atol=tol)
                @test isapprox(norm(c + A'y), 0.0, atol=tol)
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s), 0.0, atol=tol)
            end
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
        @testset "SOC3 infeasible" begin
            b = [-2, -1]
            A = [0 -1; -1 0]
            c = [0,0]
            constr_cones = [(:NonNeg,1:1),(:NonPos,2:2)]
            var_cones = [(:SOC,1:2)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Infeasible

            if duals
                y = MathProgBase.getdual(m)
                @test y[1] > 0
                @test y[2] < 0
                @test (A'y)[1]/abs((A'y)[2]) ≥ 1-tol
                @test -dot(b,y) > 0
            end
        end

        if duals
            # Problem SOC4
            # Like SOCINT1 but with copies of variables and integrality relaxed
            # Tests out-of-order indices in cones
            @testset "SOC4" begin
                b = [1.0,0.0,0.0]
                A = [1.0 0.0 0.0 0.0 0.0
                     0.0 1.0 0.0 -1.0 0.0
                     0.0 0.0 1.0 0.0 -1.0]
                c = [0.0,-2.0,-1.0,0.0,0.0]
                constr_cones = [(:Zero,[1]),(:Zero,[2,3])]
                var_cones = [(:SOC,[1,4,5]),(:Free,[2,3])]

                m = MathProgBase.ConicModel(solver)
                MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
                MathProgBase.optimize!(m)
                @test MathProgBase.status(m) == :Optimal

                x = MathProgBase.getsolution(m)
                y = MathProgBase.getdual(m)
                s = MathProgBase.getvardual(m)

                @test isapprox(dot(c,x), -dot(y,b), atol=tol)
                @test x[1]^2 ≥ x[4]^2 + x[5]^2 - tol
                @test s[1]^2 ≥ s[4]^2 + s[5]^2 - tol
                @test isapprox(s[2], 0.0, atol=tol)
                @test isapprox(s[3], 0.0, atol=tol)
                @test isapprox(norm((c+A'y) - s), 0.0, atol=tol)
            end
        end
    end
end

function conicSOCRotatedtest(solver::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    @testset "Testing SOCRotated problems through conic interface with $solver" begin
        # Problem SOCRotated1
        # min 0a + 0b - 1x - 1y
        #  st  a            == 1/2
        #  st  b            == 1
        #      2a*b >= y^2+z^2
        @testset "SOCRotated1" begin
            c = [ 0.0, 0.0, -1.0, -1.0]
            A = [ 1.0  0.0   0.0   0.0
                  0.0  1.0   0.0   0.0]
            b = [ 0.5, 1.0]
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,c,A,b,[(:Zero,1:2)],[(:SOCRotated,1:4)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 0.5, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0/sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[4], 1.0/sqrt(2.0), atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                dualobj = -dot(b,d)
                @test isapprox(dualobj, -sqrt(2.0), atol=tol)
                @test all(d .>= 0)
                vardual = MathProgBase.getvardual(m)
                @test isapprox(norm(vardual - (c+ A'd)), 0.0, atol=tol)
                @test 2*vardual[1]*vardual[2] ≥ vardual[3]^2 + vardual[4]^2 - tol
            end
        end

        # Problem SOCRotated1A - Problem SOCRotated1 with a and b substituted
        # min          -y - z
        #  st [0.5] - [      ] SOCRotated
        #     [1.0] - [      ] SOCRotated
        #     [0.0] - [-y    ] SOCRotated
        #     [0.0] - [    -z] SOCRotated
        @testset "SOCRotated1A" begin
            c = [-1.0,-1.0]
            A = [0.0 0.0; 0.0 0.0; -1.0 0.0; 0.0 -1.0]
            b = [0.5, 1.0, 0.0, 0.0]
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,c,A,b,
            [(:SOCRotated,1:4)],
            [(:Free,1:2)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0/sqrt(2.0), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 1.0/sqrt(2.0), atol=tol)
            if duals
                d = MathProgBase.getdual(m)
                @test 2*d[1]*d[2] ≥ d[3]^2 + d[4]^2 - tol
                dualobj = -dot(b,d)
                @test isapprox(dualobj, -sqrt(2.0), atol=tol)

                vardual = MathProgBase.getvardual(m)
                @test isapprox(norm(vardual - (c+ A'd)), 0.0, atol=tol)
                @test isapprox(norm(vardual), 0.0, atol=tol)
            end
        end

        # Problem SOCRotated2 - Infeasible
        # min 0
        # s.t. z ≥ 2
        #      x ≤ 1
        #      y = 1/2
        #      z^2 ≤ 2x*y
        # in conic form:
        # min 0
        # s.t. -2 + z ∈ R₊
        #      -1 + x ∈ R₋
        #     1/2 - y ∈ {0}
        #       (x,y,z) ∈ SOCRoated
        @testset "SOCRotated2 infeasible" begin
            b = [-2, -1, 1/2]
            A = [0 0 -1; -1 0 0; 0 1 0]
            c = [0,0,0]
            constr_cones = [(:NonNeg,1:1),(:NonPos,2:2),(:Zero,3:3)]
            var_cones = [(:SOCRotated,1:3)]

            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m, c, A, b, constr_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Infeasible

            if duals
                y = MathProgBase.getdual(m)
                @test y[1] > 0
                @test y[2] < 0
                vardual = MathProgBase.getvardual(m)
                @test isapprox(norm(A'y-vardual), 0.0, atol=tol)
                @test 2*vardual[1]*vardual[2] ≥ vardual[3]^2
                @test -dot(b,y) > 0
            end
        end
    end
end

function conicSOCINTtest(solver::MathProgBase.AbstractMathProgSolver;tol=1e-6)
    @testset "Testing SOCINT problems through conic interface with $solver" begin
        # Problem SOCINT1
        # min 0x - 2y - 1z
        #  st  x            == 1
        #      x >= ||(y,z)||
        #      (y,z) binary
        m = MathProgBase.ConicModel(solver)
        MathProgBase.loadproblem!(m,
        [ 0.0, -2.0, -1.0],
        [ 1.0   0.0   0.0],
        [ 1.0],
        [(:Zero,1)],
        [(:SOC,1:3)])
        MathProgBase.setvartype!(m, [:Cont,:Bin,:Bin])
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Optimal
        @test isapprox(MathProgBase.getobjval(m), -2, atol=tol)
        @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
        @test isapprox(MathProgBase.getsolution(m)[2], 1.0, atol=tol)
        @test isapprox(MathProgBase.getsolution(m)[3], 0.0, atol=tol)
    end
end


function conicEXPtest(solver::MathProgBase.AbstractMathProgSolver;duals=false, tol=1e-6)
    @testset "Testing EXP problems through conic interface with $solver" begin
        # Problem EXP1 - ExpPrimal
        # min x + y + z
        #  st  y e^(x/y) <= z, y > 0 (i.e (x, y, z) are in the exponential primal cone)
        #      x == 1
        #      y == 2
        @testset "EXP1" begin
            m = MathProgBase.ConicModel(solver)
            c = [1.0, 1.0, 1.0]
            A = [0.0 1.0 0.0;
                 1.0 0.0 0.0]
            b = [2.0, 1.0]
            MathProgBase.loadproblem!(m, c, A, b, [(:Zero,1:2)], [(:ExpPrimal, 1:3)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), (2*exp(1/2)+3), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 2.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 2*exp(1/2), atol=tol)

            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(-dot(b,d), MathProgBase.getobjval(m), atol=tol)
                u,v,w = c+A'd
                # should belong to the ExpDual cone
                @test u < 0
                @test -u*exp(v/w) ≤ exp(1)*w + tol
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s - (c+A'd)), 0.0, atol=tol)
            end

            # Permute indices - used to cause Mosek failure
            m = MathProgBase.ConicModel(solver)
            c = [1.0, 1.0, 1.0]
            A = [0.0 1.0 0.0;
                 0.0 0.0 1.0]
            b = [2.0, 1.0]
            MathProgBase.loadproblem!(m, c, A, b, [(:Zero,1:2)], [(:ExpPrimal, [3,2,1])])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), (2*exp(1/2)+3), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 2.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 2*exp(1/2), atol=tol)

            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(-dot(b,d), MathProgBase.getobjval(m), atol=tol)
                w,v,u = c+A'd
                # should belong to the ExpDual cone
                @test u < 0
                @test -u*exp(v/w) ≤ exp(1)*w + tol
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s - (c+A'd)), 0.0, atol=tol)
            end
        end

        # Problem EXP1A - ExpPrimal
        # Same as previous, except we put :ExpPrimal on constr_cones
        @testset "EXP1A" begin
            m = MathProgBase.ConicModel(solver)
            c = [1.0, 1.0, 1.0]
            A = [0.0 1.0 0.0;
                 1.0 0.0 0.0;
                -1.0 0.0 0.0;
                 0.0 -1.0 0.0;
                 0.0 0.0 -1.0]
            b = [2.0, 1.0, 0.0, 0.0, 0.0]
            MathProgBase.loadproblem!(m, c, A, b,
            [(:Zero,1:2), (:ExpPrimal, 3:5)],
            [(:Free, [2,3,1])])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), (2*exp(1/2)+3), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 2.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 2*exp(1/2), atol=tol)

            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(-dot(b,d), MathProgBase.getobjval(m), atol=tol)
                @test isapprox(norm(c+A'd), 0.0, atol=tol)
                u,v,w = d[3:5]
                # should belong to the ExpDual cone
                @test u < 0
                @test -u*exp(v/w) ≤ exp(1)*w + tol
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s), 0.0, atol=tol)
            end

            # Permute indices - used to cause Mosek failure
            m = MathProgBase.ConicModel(solver)
            c = [1.0, 1.0, 1.0]
            A = [0.0 1.0 0.0;
                -1.0 0.0 0.0;
                 1.0 0.0 0.0;
                 0.0 -1.0 0.0;
                 0.0 0.0 -1.0]
            b = [2.0, 0.0, 1.0, 0.0, 0.0]
            MathProgBase.loadproblem!(m, c, A, b,
            [(:Zero,[1,3]), (:ExpPrimal, [2,4,5])],
            [(:Free, [3,2,1])])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), (2*exp(1/2)+3), atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[1], 1.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[2], 2.0, atol=tol)
            @test isapprox(MathProgBase.getsolution(m)[3], 2*exp(1/2), atol=tol)

            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(-dot(b,d), MathProgBase.getobjval(m), atol=tol)
                @test isapprox(norm(c+A'd), 0.0, atol=tol)
                u,v,w = d[2],d[4],d[5]
                # should belong to the ExpDual cone
                @test u < 0
                @test -u*exp(v/w) ≤ exp(1)*w + tol
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s), 0.0, atol=tol)
            end
        end

        # Problem EXP2
        # A problem where ECOS was failing
        if duals
            @testset "EXP2" begin
                c = [0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
                b = [0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0]

                I = [8,11,1,4,9,12,1,4,10,13,3,7,6,7,7,8,11,14,15,9,12,14,10,13,14]
                J = [1,1,2,2,2,2,3,3,3,3,4,4,5,5,6,7,7,7,7,8,8,8,9,9,9]
                V = [-1.0,1.0,-1.0,-1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-0.5,-1.0,-0.5,1.0,-0.3,-0.3,1.0,-1.0,-0.3,-0.3,1.0,-0.3,-0.3,1.0]
                A = sparse(I, J, V, length(b), length(c))

                cone_con = [(:ExpPrimal,1:3),(:ExpPrimal,4:6),(:Zero,7:7),(:NonNeg,8:10),(:NonNeg,11:13),(:NonNeg,14:14),(:Zero,15:15)]
                cone_var = [(:Free,1:9)]
                m = MathProgBase.ConicModel(solver)
                MathProgBase.loadproblem!(m, c, A, b, cone_con, cone_var)
                MathProgBase.optimize!(m)

                primalobj = MathProgBase.getobjval(m)
                primalsol = MathProgBase.getsolution(m)
                # check primal feasibility
                conval = b-A*primalsol
                @test conval[3] + tol ≥ conval[2]*exp(conval[1]/conval[2])
                @test conval[6] + tol ≥ conval[5]*exp(conval[4]/conval[5])
                @test isapprox(conval[7], 0.0, atol=tol)
                for i in 8:14
                    @test conval[i] ≥ -tol
                end

                @test isapprox(conval[15], 0.0, atol=tol)

                y = MathProgBase.getdual(m)
                dualobj = -dot(b,y)
                @test isapprox(primalobj, dualobj, atol=tol)
                dualconval = c + A'*y
                var = MathProgBase.getvardual(m)
                @test isapprox(dualconval, var, atol=tol)
                for i in 1:9
                    @test isapprox(dualconval[i], 0.0, atol=tol)
                end
                if y[1] == 0.0
                    @test y[2] ≥ -tol
                    @test y[3] ≥ -tol
                else
                    @test y[1] < 0.0
                    @test y[3] > 0.0
                    @test -y[1]*log(-y[1]/y[3]) + y[1] - y[2] ≤ tol
                end
                if y[4] == 0.0
                    @test y[5] ≥ -tol
                    @test y[6] ≥ -tol
                else
                    @test y[4] < 0.0
                    @test y[6] > 0.0
                    @test -y[4]*log(-y[4]/y[6]) + y[4] - y[5] ≤ tol
                end
                for i in 8:14
                    @test y[i] ≥ -tol
                end
            end
        end

        # Problem EXP3
        # Another problem where ECOS was failing
        if duals
            @testset "EXP3" begin
                b = [4.0,0.0,1.0,0.0,5.0]
                c = [-1.0,0.0]

                I = [1,2,4,5]
                J = [1,1,2,2]
                V = [2.0,-1.0,-1.0,1.0]
                A = sparse(I,J,V,5,2)
                cone_var = [(:Free,[1,2])]
                cone_con = [(:NonNeg,[1]),(:ExpPrimal,[2,3,4]),(:NonNeg,[5])]
                m = MathProgBase.ConicModel(solver)
                MathProgBase.loadproblem!(m, c, A, b, cone_con, cone_var)
                MathProgBase.optimize!(m)

                primalobj = MathProgBase.getobjval(m)
                primalsol = MathProgBase.getsolution(m)
                # check primal feasibility
                conval = b-A*primalsol
                @test conval[4] + tol ≥ conval[3]*exp(conval[2]/conval[3])
                @test conval[1] ≥ -tol
                @test conval[5] ≥ -tol

                y = MathProgBase.getdual(m)
                dualobj = -dot(b,y)
                @test isapprox(primalobj, dualobj, atol=tol)
                dualconval = c + A'*y
                var = MathProgBase.getvardual(m)
                @test isapprox(norm(dualconval - var), 0.0, atol=tol)
                for i in 1:2
                    @test isapprox(dualconval[1], 0.0, atol=tol)
                end
                @test y[1] ≥ -tol
                @test -y[2]*log(-y[2]/y[4]) + y[2] - y[3] ≤ tol
                @test y[5] ≥ -tol
            end
        end

        # Problem EXPSOC1
        # A problem where MOSEK was failing
        @testset "EXPSOC1" begin
            m = MathProgBase.ConicModel(solver)
            c = [0.0, 0.0, 1.0, 0.0, 0.0]
            A = [-1.0 0.0 0.0 0.0 0.0; 1.0 3.0 1.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 2.0 3.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 1.0 0.0; -2.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 1.0 0.0 0.0 1.0]
            b = [0.0, 0.0, -1.0, 2.0, 30.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 7.0]
            con_cones = Tuple{Symbol,Array{Int64,1}}[(:NonNeg, [1]), (:Zero, [2]), (:NonNeg, [3]), (:NonNeg, [4]), (:NonNeg, [5]), (:SOC, [6, 7, 8]), (:NonNeg, [9]), (:ExpPrimal, [12, 11, 10]), (:NonNeg, [13])]
            var_cones = Tuple{Symbol,Array{Int64,1}}[(:Free, [1, 2, 3, 4, 5])]
            MathProgBase.loadproblem!(m, c, A, b, con_cones, var_cones)
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -18.082226, atol=tol)

            if duals
                d = MathProgBase.getdual(m)
                @test isapprox(-dot(b,d), MathProgBase.getobjval(m), atol=tol)
                u,v,w = d[12],d[11],d[10]
                # should belong to the ExpDual cone
                @test u < 0
                @test -u*exp(v/w) ≤ exp(1)*w + tol
            end
        end
    end
end


function conicSDPtest(solver::MathProgBase.AbstractMathProgSolver;duals=true, tol=1e-6)
    s2 = sqrt(2)
    @testset "Testing SDP problems through conic interface with $solver" begin
        # Problem SDP1 - sdo1 from MOSEK docs
        # From Mosek.jl/test/mathprogtestextra.jl, under license:
        #   Copyright (c) 2013 Ulf Worsoe, Mosek ApS
        #   Permission is hereby granted, free of charge, to any person obtaining a copy of this
        #   software and associated documentation files (the "Software"), to deal in the Software
        #   without restriction, including without limitation the rights to use, copy, modify, merge,
        #   publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
        #   to whom the Software is furnished to do so, subject to the following conditions:
        #   The above copyright notice and this permission notice shall be included in all copies or
        #   substantial portions of the Software.
        #
        #     | 2 1 0 |
        # min | 1 2 1 | . X + x1
        #     | 0 1 2 |
        #
        #
        # s.t. | 1 0 0 |
        #      | 0 1 0 | . X + x1 = 1
        #      | 0 0 1 |
        #
        #      | 1 1 1 |
        #      | 1 1 1 | . X + x2 + x3 = 1/2
        #      | 1 1 1 |
        #
        #      (x1,x2,x3) in C^3_q
        #      X in C_sdp
        #
        @testset "SDP1" begin
            m = MathProgBase.ConicModel(solver)
            #     x1   x2   x3    X11  X21  X31  X22  X32  X33
            c = [ 1.0, 0.0, 0.0,  2.0, s2,  0.0, 2.0, s2,  2.0 ]
            A = [ 1.0  0.0  0.0   1.0  0.0  0.0  1.0  0.0  1.0 ;  # A1
                  0.0  1.0  1.0   1.0  s2   s2   1.0  s2   1.0 ]  # A2
            b = [ 1.0, 0.5 ]

            MathProgBase.loadproblem!(m, c, A, b, [(:Zero,1:2)], [(:SOC,1:3),(:SDP,4:9)] )
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            pobj = MathProgBase.getobjval(m)
            @test isapprox(pobj, 7.05710509e-01, atol=tol)

            xx = MathProgBase.getsolution(m)
            x123 = xx[1:3]
            X = xx[4:9]

            if duals
                y = MathProgBase.getdual(m)
                # Check primal objective
                comp_pobj = dot(X,[2.0,s2,0.0, 2.0,s2, 2.0]) + x123[1]
                # Check dual objective
                comp_dobj = -dot(y,b)
                @test isapprox(comp_pobj, comp_dobj, atol=tol)

                var = c + A' * y
                @test (var[2]^2 + var[3]^2 - var[1]^2) < tol # (var[1],var[2],var[3]) in SOC
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(var - s), 0.0, atol=tol)

                M = [s[4]    s[5]/s2 s[6]/s2
                     s[5]/s2 s[7]    s[8]/s2
                     s[6]/s2 s[8]/s2 s[9]]

                @test eigmin(M) > -tol
            end

            # Reverse row and column indices to check nonconsecutive conic indices
            m = MathProgBase.ConicModel(solver)
            row = collect(2:-1:1)
            col = collect(9:-1:1)
            MathProgBase.loadproblem!(m, c[col], A[row,col], b[row], [(:Zero,row)], [(:SOC,col[1:3]),(:SDP,col[4:9])] )
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            pobj = MathProgBase.getobjval(m)
            @test isapprox(pobj, 7.05710509e-01, atol=tol)

            xx = MathProgBase.getsolution(m)
            x123 = xx[col[1:3]]
            X = xx[col[4:9]]

            if duals
                y = MathProgBase.getdual(m)
                # Check primal objective
                comp_pobj = dot(X,[2.0,s2,0.0, 2.0,s2, 2.0]) + x123[1]
                # Check dual objective
                comp_dobj = -dot(y,b[row])
                @test isapprox(comp_pobj, comp_dobj, atol=tol)

                var = c[col] + A[row,col]' * y
                @test (var[col[2]]^2 + var[col[3]]^2 - var[col[1]]^2) < tol # (var[1],var[2],var[3]) in SOC
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(var - s), 0.0, atol=tol)

                M = [s[col[4]]    s[col[5]]/s2 s[col[6]]/s2
                     s[col[5]]/s2 s[col[7]]    s[col[8]]/s2
                     s[col[6]]/s2 s[col[8]]/s2 s[col[9]]]

                @test eigmin(M) > -tol
            end
        end

        # Problem SDP2
        # Caused getdual to fail on SCS and Mosek
        @testset "SDP2" begin
            m = MathProgBase.ConicModel(solver)

            c = [-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-1.0]
            b = [10.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            I = [1,2,8,9,10,11,1,3,8,9,10,11,1,4,8,9,10,11,1,5,8,1,6,8,9,10,11,1,7,8,9,10,11,8,10]
            J = [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7]
            V = [1.0,1.0,-0.44999999999999996,0.45000000000000007,-0.4500000000000001,0.0,1.0,1.0,-0.7681980515339464,0.31819805153394654,-0.13180194846605373,0.0,1.0,1.0,-0.9000000000000001,0.0,0.0,0.0,1.0,1.0,-0.22500000000000003,1.0,1.0,-0.11250000000000003,0.1125,-0.11249999999999999,0.0,1.0,1.0,0.0,0.0,-0.22500000000000003,0.0,1.0,1.0]
            A = sparse(I, J, V, length(b), length(c))
            cone_con = [(:NonNeg,[1]),(:NonPos,[2,3,4,5,6,7]),(:SDP,8:10),(:Zero,11:11)]
            cone_var = [(:NonNeg,[1,2,3,4,5,6]),(:Free,[7])]

            MathProgBase.loadproblem!(m, c, A, b, cone_con, cone_var)
            MathProgBase.optimize!(m)

            @test MathProgBase.status(m) == :Optimal
            pobj = MathProgBase.getobjval(m)

            x = MathProgBase.getsolution(m)
            @test all(x[1:6] .> -tol)
            con = b - A * x
            @test eigmin([con[8] con[9]/s2 ; con[9]/s2 con[10]]) > -tol
            @test con[1] >= -tol
            @test all(con[2:7] .<= tol)
            @test isapprox(con[11], 0.0, atol=tol)

            if duals
                y = MathProgBase.getdual(m)
                @test eigmin([y[8] y[9]/s2 ; y[9]/s2 y[10]]) > -tol
                @test y[1] >= -tol
                @test all(y[2:7] .<= tol)
                var = c + A' * y
                @test all(var[1:6] .>= -tol)
                @test isapprox(var[7], 0.0, atol=tol)
                dobj = -dot(y,b)
                @test isapprox(pobj, dobj, atol=tol)
                s = MathProgBase.getvardual(m)
                @test isapprox(norm(s - var), 0.0, atol=tol)
            end
        end
    end
end
