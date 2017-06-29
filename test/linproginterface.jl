using Base.Test
using MathProgBase
using MathProgBase.SolverInterface

const MPB = MathProgBase

function test_result(m, v, c, termstatus, primalstatus, dualstatus, objval, varprim, vupdual, vlodual, conprim, condual)
    @test MPB.cangetattribute(m, MPB.TerminationStatus())
    @test MPB.getattribute(m, MPB.TerminationStatus()) == termstatus

    @test MPB.cangetattribute(m, MPB.PrimalStatus())
    @test MPB.getattribute(m, MPB.PrimalStatus()) == primalstatus

    @test MPB.cangetattribute(m, MPB.DualStatus())
    # maybe solver doesn't return dual?
    @test MPB.getattribute(m, MPB.DualStatus()) == dualstatus

    @test MPB.cangetattribute(m, MPB.ObjectiveValue())
    @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), objval, atol=eps)

    @test MPB.cangetattribute(m, MPB.VariablePrimal(), v)
    @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v), varprim, atol=eps)

    @test MPB.cangetattribute(m, MPB.VariableUpperBoundDual(), v)
    @test isapprox(MPB.getattribute(m, MPB.VariableUpperBoundDual(), v), vupdual, atol=eps)

    @test MPB.cangetattribute(m, MPB.VariableLowerBoundDual(), v)
    @test isapprox(MPB.getattribute(m, MPB.VariableLowerBoundDual(), v), vlodual, atol=eps)

    @test MPB.cangetattribute(m, MPB.ConstraintPrimal(), c)
    @test isapprox(MPB.getattribute(m, MPB.ConstraintPrimal(), c), conprim, atol=eps)

    @test MPB.cangetattribute(m, MPB.ConstraintDual(), c)
    @test isapprox(MPB.getattribute(m, MPB.ConstraintDual(), c), condual, atol=eps)
end

function linprogsolvertest(solver::AbstractMathProgSolver, eps = Base.rtoldefault(Float64))
    @testset "Testing LP interface with $solver" begin
        @testset "Basic interface" begin
            # Min -x                  Min -x
            # s.t. x + y <= 1   or    s.t. x + y - 1 ∈ NonPositive
            # x, y >= 0               x, y >= 0
            m = MPB.Model(solver)

            v = MPB.addvariables!(m, 2)

            @test MPB.cansetattribute(m, MPB.VariableLowerBound())
            MPB.setattribute!(m, MPB.VariableLowerBound(), v, [0, 0])

            @test MPB.cangetattribute(m, MPB.VariableLowerBound())
            @test isapprox(MPB.getattribute(m, MPB.VariableLowerBound(), v), [0.0, 0.0], atol=eps)

            @test MPB.cangetattribute(m, MPB.VariableUpperBound())
            @test all(MPB.getattribute(m, MPB.VariableUpperBound(), v) .>= 1e20)

            @test MPB.getattribute(m, SupportsAffineConstraint{MPB.NonPositive}())
            c = MPB.addconstraint!(m, -1, v, [1, 1], MPB.NonPositive(1))

            @test MPB.getattribute(m, MPB.VariableCount()) == 2
            @test MPB.getattribute(m, MPB.ConstraintCount()) == 1

            MPB.setattribute!(m, MPB.Sense(), MPB.MinSense)

            MPB.setobjective!(m, 0.0, v, [-1.0, 0.0])

            # TODO: query objective
            # (b, a_varref, a_coef, qi, qj, qc) = MPB.getobjective(m)
            # @test isapprox(b, 0.0, atol=eps)
            # @test a_varref == v
            # @test isapprox(a_coef, [-1.0, 0.0], atol=eps)

            MPB.optimize!(m)

            test_result(m, v, c, MPB.Success, MPB.FeasiblePoint, MPB.FeasiblePoint, -1, [1.0, 0.0], [0.0, 0.0], [0.0, 1.0], [1.0], [-1.0])

            # change objective to Max +x
            @test MPB.cansetattribute(m, MPB.Sense())
            MPB.setattribute!(m, MPB.Sense(), MPB.MaxSense)
            MPB.modifyobjective!(m, 1, v[1], 1.0)

            MPB.optimize!(m)

            test_result(m, v, c, MPB.Success, MPB.FeasiblePoint, MPB.FeasiblePoint, 1, [1.0, 0.0], [0.0, -1.0], [0.0, 0.0], [1.0], [1.0])

        end
    end
end
"""
            # add new variable to get:
            # max x + 2z
            # s.t. x + y + z <= 1
            # x,y,z >= 0
            addvar!(m, [1], [1.0], 0, Inf, 2.0)

            @test numvar(m) == 3
            @test numconstr(m) == 1

            optimize!(m)

            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 2, atol=eps)
            @test isapprox(norm(getsolution(m) - [0.0, 0.0, 1.0]), 0.0, atol=eps)
            @test isapprox(getconstrsolution(m)[1], 1.0, atol=eps)
            @test isapprox(getconstrduals(m)[1], 2.0, atol=eps)
            @test isapprox(norm(getreducedcosts(m) - [-1.0, -2.0, 0.0]), 0.0, atol=eps)

            setvarLB!(m, [-1.0,0.0,0.0])
            optimize!(m)

            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 3, atol=eps)

            # fix z to zero
            setvarLB!(m, [0.0,0.0,0.0])
            setvarUB!(m, [Inf,Inf,0.0])
            optimize!(m)

            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 1, atol=eps)

            setconstrUB!(m, [2.0])
            setconstrLB!(m, [2.0])
            optimize!(m)
            @test isapprox(getobjval(m), 2, atol=eps)

            setobj!(m, [1.0, 2.0, 0.0])
            optimize!(m)
            @test isapprox(getobjval(m), 4, atol=eps)
            @test isapprox(norm(getsolution(m) - [0.0, 2.0, 0.0]), 0.0, atol=eps)


            # we now have:
            # max x+2y
            # s.t. x + y + z == 2
            # x,y >= 0, z = 0

            # add constraint x - y >= 0
            addconstr!(m, [1,2], [1.0,-1.0], 0.0, Inf)
            optimize!(m)

            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 3, atol=eps)
            @test isapprox(norm(getsolution(m) - [1.0, 1.0, 0.0]), 0.0, atol=eps)
            @test isapprox(norm(getconstrduals(m) - [1.5, -0.5]), 0.0, atol=eps)
            @test isapprox(norm(getreducedcosts(m) - [0.0, 0.0, -1.5]), 0.0, atol=eps)
        end

"""

function linprogsolvertest2(solver::AbstractMathProgSolver, eps = Base.rtoldefault(Float64))
    @testset "Testing LP interface with $solver part 2" begin
        @testset "addvar! interface" begin
            m = LinearQuadraticModel(solver)

            # Min -x
            # s.t. x + y <= 1
            # x, y >= 0
            m = MPB.Model(solver)

            x = MPB.addvariable!(m)
            y = MPB.addvariable!(m)

            MPB.setattribute!(m, MPB.VariableLowerBound(), [x, y], [0, 0])

            c = MPB.addconstraint!(m, -1, [x, y], [1.0, 1.0], MPB.NonPositive(1))

            MPB.setattribute!(m, MPB.Sense(), MPB.MinSense)
            MPB.setobjective!(m, 0.0, x, -1.0)

            @test MPB.getattribute(m, MPB.VariableCount()) == 2
            @test MPB.getattribute(m, MPB.ConstraintCount()) == 1

            # getters #190
            # @test getvarLB(m) == [0.,0.]
            # @test getvarUB(m)[1] >= 1e20
            # @test getvarUB(m)[2] >= 1e20
            # @test getconstrLB(m)[1] <= -1e20
            # @test getconstrUB(m)[1] == 1.0
            # @test getobj(m) == [-1.,0.]
            # @test getsense(m) == :Min

            optimize!(m)
            @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
            @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
            @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

            @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), -1.0, atol=eps)
            @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), [x, y]), [1.0, 0.0], atol=eps)

            @test isapprox(MPB.getattribute(m, MPB.ConstraintPrimal()), 0.0, atol=eps)

            @test isapprox(MPB.getattribute(m, MPB.ConstraintDual()), -1.0, atol=eps)
            @test isapprox(MPB.getattribute(m, MPB.VariableDual(), [x, y]), [0.0, 1.0], atol=eps)
        end

        @testset "setconstrLB! and setconstrUB!" begin
            ####################################
            # test setconstrLB! and setconstrUB!
            # Test that:
            #   Modifying lower bound works
            #   Modifying upper bound works
            #   Setting upper and lower bound to same value works
            #

            m = MPB.Model(solver)

            x = MPB.addvariable!(m)
            y = MPB.addvariable!(m)

            # Min  x - y
            # s.t. 0.0 <= x          (c1)
            #             y <= 0.0   (c2)
            # x,y unbounded
            c1 = MPB.addconstraint!(m, 0.0, x, [1.0, 0.0], MPB.NonNegative(1))
            c2 = MPB.addconstraint!(m, 0.0, y, [0.0, 1.0], MPB.NonPositive(1))

            MPB.setattribute!(m, MPB.Sense(), MPB.MinSense)

            MPB.setobjective!(m, 0.0, [x, y], [1.0, -1.0])

            MPB.optimize!(m)

            @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
            @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
            @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint
            @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 0.0, atol=eps)
            @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), [x, y]), [0.0, 0.0], atol=eps)

            # Min  x - y
            # s.t. 100.0 <= x          or -100 + x \in NonNegative
            #               y <= 0.0
            # x,y unbounded
            modifyconstraint!(m, c1, 1, -100.0)
            optimize!(m)
            @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
            @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
            @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint
            @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 100.0, atol=eps)
            @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), [x, y]), [100.0, 0.0], atol=eps)

            # Min  x - y
            # s.t. 100.0 <= x
            #               y <= -100.0
            # x,y unbounded
            modifyconstraint!(m, c2, 1, 100.0)
            optimize!(m)
            @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
            @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
            @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint
            @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 200.0, atol=eps)
            @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), [x, y]), [100.0, -100.0], atol=eps)
        end

        @testset "Issue #40 from Gurobi.jl" begin
            # Test issue #40 from Gurobi.jl
            # min  x
            # s.t. x >= 0 (c1)
            #      x >= 3 (c2)

            m = MPB.Model(solver)

            x = MPB.addvariable!(m)

            c = MPB.addconstraint!(m, [0.0 -3.0], [x x], [1.0 1.0], MPB.NonNegative(2))

            MPB.setattribute!(m, MPB.Sense(), MPB.MinSense)

            MPB.setobjective!(m, 0.0, [x], [1.0])

            MPB.optimize!(m)

            for i = 1:MPB.getattribute(m, MPB.ConstraintCount())
                # Get Bounds #190
                # @test getconstrLB(m)[i] <= getconstrsolution(m)[i] + eps
                # @test getconstrsolution(m)[i] <= getconstrUB(m)[i] + eps
            end

            # min  x
            # s.t. x <= 0
            #      x <= 3
            m = MPB.Model(solver)

            x = MPB.addvariable!(m)

            c = MPB.addconstraint!(m, [0.0 -3.0], [x x], [1.0 1.0], MPB.NonPositive(2))

            MPB.setattribute!(m, MPB.Sense(), MPB.MinSense)

            MPB.setobjective!(m, 0.0, [x], [1.0])

            MPB.optimize!(m)
            for i = 1:MPB.getattribute(m, MPB.ConstraintCount())
                # Get Bounds #190
                # @test getconstrLB(m)[i] <= getconstrsolution(m)[i] + eps
                # @test getconstrsolution(m)[i] <= getconstrUB(m)[i] + eps
            end
        end

        @testset "Change coeffs, del constr, del var" begin
            #####################################
            # Start from simple LP
            # Solve it
            # Copy and solve again
            # Chg coeff, solve, change back solve
            # del constr and solve
            # del var and solve

            #   maximize x + y
            #
            #   s.t. 2 x + 1 y <= 4 (c1)
            #        1 x + 2 y <= 4 (c2)
            #        x >= 0, y >= 0
            #
            #   solution: x = 1.3333333, y = 1.3333333, objv = 2.66666666

            m = MPB.Model(solver)

            v = MPB.addvariables!(m, 2)
            
            MPB.setattribute!(m, MPB.VariableLowerBound(), v, [0, 0])

            c1 = MPB.addconstraint!(m, -4, v, [2, 1], MPB.NonPositive(1))
            c2 = MPB.addconstraint!(m, -4, v, [1, 2], MPB.NonPositive(1))

            MPB.setattribute!(m, MPB.Sense(), MPB.MaxSense)

            MPB.setobjective!(m, 0.0, v, [1.0, 1.0])

            MPB.optimize!(m)

            @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
            @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
            @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

            @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 2.6666666666, atol=eps)

            @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v), [1.3333333333, 1.3333333333], atol=eps)

            # copy and solve again
            
            if applicable(copy, m) # see 188
                m2 = copy(m)
                
                optimize!(m2)

                @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
                @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
                @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

                @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 2.6666666666, atol=eps)
                @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v), [1.3333333333, 1.3333333333], atol=eps)
            end


            # change coeff
            #   maximize x + y
            #
            #   s.t. 2 x + 2 y <= 4
            #        1 x + 2 y <= 4
            #        x >= 0, y >= 0
            #
            #   solution: x = 0, y = 2, objv = 2
            if applicable(changecoeffs!, m, [1], [2],  [2.])
                modifyconstraint!(m, c1, v[2], 2.0)
                optimize!(m)

                MPB.optimize!(m)

                @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
                @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
                @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

                @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 2.0, atol=eps)

                @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v), [0.0, 2.0], atol=eps)
            end


            # delconstrs and solve
            #   maximize x + y
            #
            #   s.t. 1 x + 2 y <= 4
            #        x >= 0, y >= 0
            #
            #   solution: x = 4, y = 0, objv = 4
            if applicable(delete!, m, [1])
                delete!(m, c1)

                optimize!(m)
                
                @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
                @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
                @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

                @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 4.0, atol=eps)

                @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v), [4.0, 0.0], atol=eps)
            end

            # delvars and solve
            #   maximize y
            #
            #   s.t.  2 y <= 4
            #           y >= 0
            #
            #   solution: y = 2, objv = 2
            if applicable(delete!, m, [1])
                delete!(m, v[1])

                optimize!(m)

                @test MPB.getattribute(m, MPB.TerminationStatus()) == :Success
                @test MPB.getattribute(m, MPB.PrimalStatus()) == MPB.FeasiblePoint
                @test MPB.getattribute(m, MPB.DualStatus()) == MPB.FeasiblePoint

                @test isapprox(MPB.getattribute(m, MPB.ObjectiveValue()), 2.0, atol=eps)

                @test isapprox(MPB.getattribute(m, MPB.VariablePrimal(), v[2]), 2.0, atol=eps)
            end
        end
    end
end

"""

function linprogsolvertestextra(solver::AbstractMathProgSolver; eps = Base.rtoldefault(Float64))
    @testset "Testing LP interface extra with $solver" begin
        ####################################
        # test setconstrLB! and setconstrUB!
        # Test that:
        #   Modifying lower bound works
        #   Modifying upper bound works
        #   Setting upper and lower bound to same value works
        #
        @testset "setconstrLB! and setcontrUB!" begin
            m = LinearQuadraticModel(solver)
            # Min  x - y
            # s.t. 0.0 <= x <= 0.0
            #      0.0 <= y <= 0.0
            # x,y unbounded
            loadproblem!(m, [ 1.0 0.0 ; 0.0 1.0 ], [-Inf, -Inf], [Inf,Inf], [1.0, -1.0], [0.0,0.0], [0.0,0.0], :Min)

            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m)), 0.0)

            # Min  x - y
            # s.t. 0.0 <= x <= 100.0
            #      0.0 <= y <= 100.0
            # x,y unbounded
            setconstrUB!(m,[100.0,100.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), -100.0)
            @test isapprox(norm(getsolution(m) - [ 0.0, 100.0 ]), 0.0)

            # Min  x - y
            # s.t. -100.0 <= x <= 100.0
            #      -100.0 <= y <= 100.0
            # x,y unbounded
            setconstrLB!(m,[-100.0,-100.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), -200.0)
            @test isapprox(norm(getsolution(m) - [ -100.0, 100.0 ]), 0.0)

            # Min  x - y
            # s.t. -100.0 <= x <= 100.0
            #      -100.0 <= y <= 100.0
            # x,y unbounded
            setconstrLB!(m,[10.0,10.0])
            setconstrUB!(m,[10.0,10.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m) - [ 10.0, 10.0 ]), 0.0)

            # Min  x - y
            # s.t. 0.0  <= x <= Inf
            #      -Inf <= y <= 0.0
            # x,y unbounded
            setconstrLB!(m,[0.0,-Inf])
            setconstrUB!(m,[Inf,0.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m)), 0.0)

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= 0.0
            # x,y unbounded
            setconstrLB!(m,[-Inf,-Inf])
            setconstrUB!(m,[Inf,0.0])
            optimize!(m)
            @test status(m) == :Unbounded
            @test dot(getunboundedray(m),[1.0,-1.0]) < 1e-8


            # Min  x - y
            # s.t. 0.0 <= x <= Inf
            #      -Inf <= y <= Inf
            # x,y unbounded
            setconstrLB!(m,[0.0,-Inf])
            setconstrUB!(m,[Inf,Inf])
            optimize!(m)
            @test status(m) == :Unbounded
            @test dot(getunboundedray(m),[1.0,-1.0]) < 1e-8
        end

        ####################################
        # test setvarLB! and setvarUB!
        # Test that:
        #   Modifying lower bound works
        #   Modifying upper bound works
        #   Setting upper and lower bound to same value works
        @testset "setvarLB! and setvarUB!" begin
            m = LinearQuadraticModel(solver)
            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # 0.0 <= x <= 0.0
            # 0.0 <= y <= 0.0
            loadproblem!(m, [ 1.0 0.0 ; 0.0 1.0 ], [0.0,0.0], [0.0,0.0], [1.0, -1.0], [-Inf, -Inf], [Inf,Inf], :Min)

            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m)), 0.0)



            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # 0.0 <= x <= 100.0
            # 0.0 <= y <= 100.0

            # x,y unbounded
            setvarUB!(m,[100.0,100.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), -100.0)
            @test isapprox(norm(getsolution(m) - [ 0.0, 100.0 ]), 0.0)

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # -100.0 <= x <= 100.0
            # -100.0 <= y <= 100.0
            setvarLB!(m,[-100.0,-100.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), -200.0)
            @test isapprox(norm(getsolution(m) - [ -100.0, 100.0 ]), 0.0)

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # -100.0 <= x <= 100.0
            # -100.0 <= y <= 100.0
            setvarLB!(m,[10.0,10.0])
            setvarUB!(m,[10.0,10.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m) - [ 10.0, 10.0 ]), 0.0)

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # 0.0  <= x <= Inf
            # -Inf <= y <= 0.0
            # x,y unbounded
            setvarLB!(m,[0.0,-Inf])
            setvarUB!(m,[Inf,0.0])
            optimize!(m)
            @test status(m) == :Optimal
            @test isapprox(getobjval(m), 0.0)
            @test isapprox(norm(getsolution(m)), 0.0)

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            # -Inf <= x <= Inf
            # -Inf <= y <= 0.0
            setvarLB!(m,[-Inf,-Inf])
            setvarUB!(m,[Inf,0.0])
            optimize!(m)
            @test status(m) == :Unbounded
            @test dot(getunboundedray(m),[1.0,-1.0]) < 1e-8

            # Min  x - y
            # s.t. -Inf <= x <= Inf
            #      -Inf <= y <= Inf
            #  0.0 <= x <= Inf
            # -Inf <= y <= Inf
            # x,y unbounded
            setvarLB!(m,[0.0,-Inf])
            setvarUB!(m,[Inf,Inf])
            optimize!(m)
            @test status(m) == :Unbounded
            @test dot(getunboundedray(m),[1.0,-1.0]) < 1e-8
        end
    end
end
"""
