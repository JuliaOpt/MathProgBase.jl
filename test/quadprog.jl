using Base.Test
using MathProgBase
using MathProgBase.SolverInterface

function quadprogtest(solver=MathProgBase.defaultQPsolver)
    println("Testing quadprog with solver ", string(typeof(solver)))

    sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf,solver)
    @test sol.status == :Optimal
    @test_approx_eq_eps sol.objval 130/70 1e-6
    @test_approx_eq_eps norm(sol.sol[1:3] - [0.5714285714285715,0.4285714285714285,0.8571428571428572]) 0.0 1e-6

    let
    m = LinearQuadraticModel(solver)
    loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)

    setquadobj!(m,diagm([10.0,10.0,10.0]))
    rows = [1, 2, 2, 2, 3, 3, 3]
    cols = [1, 1, 1, 2, 2, 3, 3]
    vals = Float64[2, 0.5, 0.5, 2, 1, 1, 1]
    setquadobj!(m,rows,cols,vals)
    optimize!(m)
    stat = status(m)
    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) 130/70 1e-6
    @test_approx_eq_eps norm(getsolution(m) - [0.5714285714285715,0.4285714285714285,0.8571428571428572]) 0.0 1e-6
    end

    let
    m = LinearQuadraticModel(solver)
    loadproblem!(m, [-1. 1.; 1. 1.], [0.,0.], [Inf,Inf], [1.,1.], [0.,0.], [Inf,Inf], :Max)
    addquadconstr!(m, [2], [1.], [1], [1], [1.], '<', 2)
    optimize!(m)
    stat = status(m)
    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) 2.25 1e-6
    @test_approx_eq_eps norm(getsolution(m) - [0.5,1.75]) 0.0 1e-3
    end
    println("Done")
end

function qpdualtest(solver=MathProgBase.defaultQPsolver)
    println("Testing QP duals with solver ", string(typeof(solver)))
    # max x
    # s.t. x^2 <= 2
    m = LinearQuadraticModel(solver)
    loadproblem!(m, Array(Float64,0,1), [-Inf], [Inf], [1.0], Float64[], Float64[], :Max)
    addquadconstr!(m, [], [], [1], [1], [1.0], '<', 2.0)
    optimize!(m)
    stat = status(m)

    @test numlinconstr(m) == 0
    @test numquadconstr(m) == 1
    @test numconstr(m) == 1
    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) sqrt(2) 1e-6
    @test_approx_eq_eps getsolution(m)[1] sqrt(2) 1e-6
    @test_approx_eq_eps getquadconstrduals(m)[1] 0.5/sqrt(2) 1e-6

    # min -x
    # s.t. x^2 <= 2
    m = LinearQuadraticModel(solver)
    loadproblem!(m, Array(Float64,0,1), [-Inf], [Inf], [-1.0], Float64[], Float64[], :Min)
    addquadconstr!(m, [], [], [1], [1], [1.0], '<', 2.0)
    optimize!(m)
    stat = status(m)

    @test numlinconstr(m) == 0
    @test numquadconstr(m) == 1
    @test numconstr(m) == 1
    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) -sqrt(2) 1e-6
    @test_approx_eq_eps getsolution(m)[1] sqrt(2) 1e-6
    @test_approx_eq_eps getquadconstrduals(m)[1] -0.5/sqrt(2) 1e-6

    println("Done")
end

function socptest(solver=MathProgBase.defaultQPsolver)
    println("Testing SOCP interface with solver ", string(typeof(solver)))
    
    # min t
    # s.t. x + y >= 1
    #      x^2 + y^2 <= t^2
    #      t >= 0
    m = LinearQuadraticModel(solver)
    loadproblem!(m, [ 1.0 1.0 0.0 ], [-Inf,-Inf,0.0], [Inf,Inf,Inf], [0.0,0.0,1.0], [1.0],[Inf], :Min)
    addquadconstr!(m, [], [], [1,2,3], [1,2,3], [1.0,1.0,-1.0],'<',0.0)
    optimize!(m)
    stat = status(m)

    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) sqrt(1/2) 1e-6
    @test_approx_eq_eps norm(getsolution(m) - [0.5,0.5,sqrt(1/2)]) 0.0 1e-3
    println("Done")


end
