using Base.Test
using MathProgBase
using MathProgBase.SolverInterface

function linprogsolvertest(solver::AbstractMathProgSolver, eps = Base.rtoldefault(Float64))

    m = LinearQuadraticModel(solver)

    # Min -x
    # s.t. x + y <= 1
    # x, y >= 0
    loadproblem!(m, [1. 1.], [0.,0.], [Inf, Inf], [-1.,0.], [-Inf],[1.], :Min)
    @test numvar(m) == 2
    @test numconstr(m) == 1
    @test getvarLB(m) == [0.,0.]
    @test getvarUB(m)[1] >= 1e20
    @test getvarUB(m)[2] >= 1e20
    @test getconstrLB(m)[1] < -1e30
    @test getconstrUB(m)[1] == 1.0
    @test getobj(m) == [-1.,0.]
    @test getsense(m) == :Min

    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) -1 eps
    @test_approx_eq_eps getsolution(m) [1.0, 0.0] eps
    @test_approx_eq_eps getconstrsolution(m) [1.0] eps
    @test_approx_eq_eps getconstrduals(m) [-1.0] eps
    @test_approx_eq_eps getreducedcosts(m) [0.0, 1.0] eps

    setsense!(m, :Max)
    # max x
    setobj!(m, [1.0,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 1 eps
    @test_approx_eq_eps getsolution(m) [1.0, 0.0] eps
    @test_approx_eq_eps getconstrsolution(m) [1.0] eps
    @test_approx_eq_eps getconstrduals(m) [1.0] eps
    @test_approx_eq_eps getreducedcosts(m) [0.0, -1.0] eps

    # add new variable to get:
    # max x + 2z
    # s.t. x + y + z <= 1
    # x,y,z >= 0
    addvar!(m, [1], [1.0], 0, Inf, 2.0)

    @test numvar(m) == 3
    @test numconstr(m) == 1

    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 2 eps
    @test_approx_eq_eps getsolution(m) [0.0, 0.0, 1.0] eps
    @test_approx_eq_eps getconstrsolution(m) [1.0] eps
    @test_approx_eq_eps getconstrduals(m) [2.0] eps
    @test_approx_eq_eps getreducedcosts(m) [-1.0, -2.0, 0.0] eps

    setvarLB!(m, [-1.0,0.0,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 3 eps

    # fix z to zero
    setvarLB!(m, [0.0,0.0,0.0])
    setvarUB!(m, [Inf,Inf,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 1 eps

    setconstrUB!(m, [2.0])
    setconstrLB!(m, [2.0])
    optimize!(m)
    @test_approx_eq_eps getobjval(m) 2 eps

    setobj!(m, [1.0, 2.0, 0.0])
    optimize!(m)
    @test_approx_eq_eps getobjval(m) 4 eps
    @test_approx_eq_eps getsolution(m) [0.0, 2.0, 0.0] eps


    # we now have:
    # max x+2y
    # s.t. x + y + z == 2
    # x,y >= 0, z = 0

    # add constraint x - y >= 0
    addconstr!(m, [1,2], [1.0,-1.0], 0.0, Inf)
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 3 eps
    @test_approx_eq_eps getsolution(m) [1.0, 1.0, 0.0] eps
    @test_approx_eq_eps getconstrduals(m) [1.5, -0.5] eps
    @test_approx_eq_eps getreducedcosts(m) [0.0, 0.0, -1.5] eps


    # test addvar! interface

    m = LinearQuadraticModel(solver)

    # Min -x
    # s.t. x + y <= 1
    # x, y >= 0
    addvar!(m, 0, Inf, -1)
    addvar!(m, 0, Inf, 0)
    addconstr!(m, [1, 2], [1.0, 1.0], -Inf, 1.0)
    @test numvar(m) == 2
    @test numconstr(m) == 1
    @test getvarLB(m) == [0.,0.]
    @test getvarUB(m)[1] >= 1e20
    @test getvarUB(m)[2] >= 1e20
    @test getconstrLB(m)[1] <= -1e20
    @test getconstrUB(m)[1] == 1.0
    @test getobj(m) == [-1.,0.]
    @test getsense(m) == :Min

    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) -1 eps
    @test_approx_eq_eps getsolution(m) [1.0, 0.0] eps
    @test_approx_eq_eps getconstrsolution(m) [1.0] eps
    @test_approx_eq_eps getconstrduals(m) [-1.0] eps
    @test_approx_eq_eps getreducedcosts(m) [0.0, 1.0] eps


    ####################################
    # test setconstrLB! and setconstrUB!
    # Test that:
    #   Modifying lower bound works
    #   Modifying upper bound works
    #   Setting upper and lower bound to same value works
    #

    m = LinearQuadraticModel(solver)
    # Min  x - y
    # s.t. 0.0 <= x
    #             y <= 0.0
    # x,y unbounded

    loadproblem!(m, [ 1.0 0.0 ; 0.0 1.0 ], [-Inf, -Inf], [Inf,Inf], [1.0, -1.0], [0.0,-Inf], [Inf,0.0], :Min)

    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 0.0 eps
    @test_approx_eq_eps getsolution(m) [ 0.0, 0.0 ] eps

    # Min  x - y
    # s.t. 100.0 <= x
    #               y <= 0.0
    # x,y unbounded
    setconstrLB!(m,[100.0,-Inf])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 100.0 eps
    @test_approx_eq_eps getsolution(m) [ 100.0, 0.0 ] eps

    # Min  x - y
    # s.t. 100.0 <= x
    #               y <= -100.0
    # x,y unbounded
    setconstrUB!(m,[Inf,-100.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 200.0 eps
    @test_approx_eq_eps getsolution(m) [ 100.0, -100.0 ] eps

    # Test issue #40 from Gurobi.jl
    # min  x
    # s.t. x >= 0
    #      x >= 3
    m = LinearQuadraticModel(solver)
    loadproblem!(m, [1.0 1.0]', [-Inf], [Inf], [1.0], [0.0, 3.0], [Inf, Inf], :Min)
    optimize!(m)
    for i = 1:length(getconstrLB(m))
        @test getconstrLB(m)[i] <= getconstrsolution(m)[i] + eps
        @test getconstrsolution(m)[i] <= getconstrUB(m)[i] + eps
    end

    # min  x
    # s.t. x <= 0
    #      x <= 3
    m = LinearQuadraticModel(solver)
    loadproblem!(m, [1.0 1.0]', [-Inf], [Inf], [1.0], [-Inf, -Inf], [0.0, 3.0], :Max)
    optimize!(m)
    for i = 1:length(getconstrLB(m))
        @test getconstrLB(m)[i] <= getconstrsolution(m)[i] + eps
        @test getconstrsolution(m)[i] <= getconstrUB(m)[i] + eps
    end


    #####################################
    # Start from simple LP
    # Solve it
    # Copy and solve again
    # Chg coeff, solve, change back solve
    # del constr and solve
    # del var and solve

    #   maximize x + y
    #
    #   s.t. 2 x + 1 y <= 4
    #        1 x + 2 y <= 4
    #        x >= 0, y >= 0
    #
    #   solution: x = 1.3333333, y = 1.3333333, objv = 2.66666666

    m = LinearQuadraticModel(solver)

    loadproblem!(m, [ 2.0 1.0 ; 1.0 2.0 ], [0.0,0.0], [Inf,Inf], [1.0, 1.0], [-Inf, -Inf], [4.0,4.0], :Max)

    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq_eps getobjval(m) 2.6666666666 eps
    @test_approx_eq_eps getsolution(m) [1.3333333333, 1.3333333333] eps

    # copy and solve again

    if applicable(copy, m)
        m2 = copy(m)
        optimize!(m2)
        @test status(m2) == :Optimal
        @test_approx_eq_eps getobjval(m2) 2.6666666666 eps
        @test_approx_eq_eps getsolution(m2) [1.3333333333, 1.3333333333] eps
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
        changecoeffs!(m, [1], [2],  [2.])
        optimize!(m)
        @test status(m) == :Optimal
        @test_approx_eq_eps getobjval(m) 2.0 eps
        @test_approx_eq_eps getsolution(m) [0.0, 2.0] eps
    end


    # delconstrs and solve
    #   maximize x + y
    #
    #   s.t. 1 x + 2 y <= 4
    #        x >= 0, y >= 0
    #
    #   solution: x = 4, y = 0, objv = 4
    if applicable(delconstrs!, m, [1])
        delconstrs!(m, [1])
        optimize!(m)
        @test status(m) == :Optimal
        @test_approx_eq_eps getobjval(m) 4.0 eps
        @test_approx_eq_eps getsolution(m) [4.0, 0.0] eps
    end

    # delvars and solve
    #   maximize y
    #
    #   s.t.  2 y <= 4
    #           y >= 0
    #
    #   solution: y = 2, objv = 2
    if applicable(delvars!, m, [1])
        delvars!(m, [1])
        optimize!(m)
        @test status(m) == :Optimal
        @test_approx_eq_eps getobjval(m) 2.0 eps
        @test_approx_eq_eps getsolution(m) [2.0] eps
    end


end



function linprogsolvertestextra(solver::AbstractMathProgSolver; eps = Base.rtoldefault(Float64))
    ####################################
    # test setconstrLB! and setconstrUB!
    # Test that:
    #   Modifying lower bound works
    #   Modifying upper bound works
    #   Setting upper and lower bound to same value works
    #

    m = LinearQuadraticModel(solver)
    # Min  x - y
    # s.t. 0.0 <= x <= 0.0
    #      0.0 <= y <= 0.0
    # x,y unbounded
    loadproblem!(m, [ 1.0 0.0 ; 0.0 1.0 ], [-Inf, -Inf], [Inf,Inf], [1.0, -1.0], [0.0,0.0], [0.0,0.0], :Min)

    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 0.0, 0.0 ]

    # Min  x - y
    # s.t. 0.0 <= x <= 100.0
    #      0.0 <= y <= 100.0
    # x,y unbounded
    setconstrUB!(m,[100.0,100.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) -100.0
    @test_approx_eq getsolution(m) [ 0.0, 100.0 ]

    # Min  x - y
    # s.t. -100.0 <= x <= 100.0
    #      -100.0 <= y <= 100.0
    # x,y unbounded
    setconstrLB!(m,[-100.0,-100.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) -200.0
    @test_approx_eq getsolution(m) [ -100.0, 100.0 ]

    # Min  x - y
    # s.t. -100.0 <= x <= 100.0
    #      -100.0 <= y <= 100.0
    # x,y unbounded
    setconstrLB!(m,[10.0,10.0])
    setconstrUB!(m,[10.0,10.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 10.0, 10.0 ]

    # Min  x - y
    # s.t. 0.0  <= x <= Inf
    #      -Inf <= y <= 0.0
    # x,y unbounded
    setconstrLB!(m,[0.0,-Inf])
    setconstrUB!(m,[Inf,0.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 0.0, 0.0 ]

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





    ####################################
    # test setvarLB! and setvarUB!
    # Test that:
    #   Modifying lower bound works
    #   Modifying upper bound works
    #   Setting upper and lower bound to same value works
    

    m = LinearQuadraticModel(solver)
    # Min  x - y
    # s.t. -Inf <= x <= Inf
    #      -Inf <= y <= Inf
    # 0.0 <= x <= 0.0
    # 0.0 <= y <= 0.0
    loadproblem!(m, [ 1.0 0.0 ; 0.0 1.0 ], [0.0,0.0], [0.0,0.0], [1.0, -1.0], [-Inf, -Inf], [Inf,Inf], :Min)

    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 0.0, 0.0 ]



    # Min  x - y
    # s.t. -Inf <= x <= Inf
    #      -Inf <= y <= Inf
    # 0.0 <= x <= 100.0
    # 0.0 <= y <= 100.0

    # x,y unbounded
    setvarUB!(m,[100.0,100.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) -100.0
    @test_approx_eq getsolution(m) [ 0.0, 100.0 ]

    # Min  x - y
    # s.t. -Inf <= x <= Inf
    #      -Inf <= y <= Inf
    # -100.0 <= x <= 100.0
    # -100.0 <= y <= 100.0
    setvarLB!(m,[-100.0,-100.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) -200.0
    @test_approx_eq getsolution(m) [ -100.0, 100.0 ]

    # Min  x - y
    # s.t. -Inf <= x <= Inf
    #      -Inf <= y <= Inf
    # -100.0 <= x <= 100.0
    # -100.0 <= y <= 100.0
    setvarLB!(m,[10.0,10.0])
    setvarUB!(m,[10.0,10.0])
    optimize!(m)
    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 10.0, 10.0 ]

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
    @test_approx_eq getobjval(m) 0.0
    @test_approx_eq getsolution(m) [ 0.0, 0.0 ]

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
