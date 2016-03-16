using Base.Test
using MathProgBase
using MathProgBase.SolverInterface

function linprogsolvertest(solver::AbstractMathProgSolver)

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
    @test_approx_eq getobjval(m) -1
    @test_approx_eq getsolution(m) [1.0, 0.0]
    @test_approx_eq getconstrsolution(m) [1.0]
    @test_approx_eq getconstrduals(m) [-1.0]
    @test_approx_eq getreducedcosts(m) [0.0, 1.0]

    setsense!(m, :Max)
    # max x
    setobj!(m, [1.0,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 1
    @test_approx_eq getsolution(m) [1.0, 0.0]
    @test_approx_eq getconstrsolution(m) [1.0]
    @test_approx_eq getconstrduals(m) [1.0]
    @test_approx_eq getreducedcosts(m) [0.0, -1.0]

    # add new variable to get:
    # max x + 2z
    # s.t. x + y + z <= 1
    # x,y,z >= 0
    addvar!(m, [1], [1.0], 0, Inf, 2.0)

    @test numvar(m) == 3
    @test numconstr(m) == 1

    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 2
    @test_approx_eq getsolution(m) [0.0, 0.0, 1.0]
    @test_approx_eq getconstrsolution(m) [1.0]
    @test_approx_eq getconstrduals(m) [2.0]
    @test_approx_eq getreducedcosts(m) [-1.0, -2.0, 0.0]

    setvarLB!(m, [-1.0,0.0,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 3

    # fix z to zero
    setvarLB!(m, [0.0,0.0,0.0])
    setvarUB!(m, [Inf,Inf,0.0])
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 1

    setconstrUB!(m, [2.0])
    setconstrLB!(m, [2.0])
    optimize!(m)
    @test_approx_eq getobjval(m) 2

    setobj!(m, [1.0, 2.0, 0.0])
    optimize!(m)
    @test_approx_eq getobjval(m) 4
    @test_approx_eq getsolution(m) [0.0, 2.0, 0.0]


    # we now have:
    # max x+2y
    # s.t. x + y + z == 2
    # x,y >= 0, z = 0

    # add constraint x - y >= 0
    addconstr!(m, [1,2], [1.0,-1.0], 0.0, Inf)
    optimize!(m)

    @test status(m) == :Optimal
    @test_approx_eq getobjval(m) 3
    @test_approx_eq getsolution(m) [1.0, 1.0, 0.0]
    @test_approx_eq getconstrduals(m) [1.5, -0.5]
    @test_approx_eq getreducedcosts(m) [0.0, 0.0, -1.5]


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
    @test_approx_eq getobjval(m) -1
    @test_approx_eq getsolution(m) [1.0, 0.0]
    @test_approx_eq getconstrsolution(m) [1.0]
    @test_approx_eq getconstrduals(m) [-1.0]
    @test_approx_eq getreducedcosts(m) [0.0, 1.0]


end
