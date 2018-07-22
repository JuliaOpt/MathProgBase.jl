using Compat.Test
using Compat.LinearAlgebra
using MathProgBase

function quadprogtest(solver)
    @testset "Testing quadprog with $solver" begin
        sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf,solver)
        @test sol.status == :Optimal
        @test isapprox(sol.objval, 130/70, atol=1e-6)
        @test isapprox(norm(sol.sol[1:3] - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)

        @testset "QP1" begin
            m = MathProgBase.LinearQuadraticModel(solver)
            MathProgBase.loadproblem!(m, [1. 2. 3.; 1. 1. 0.],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[0.,0.,0.],[4., 1.],[Inf,Inf], :Min)

            MathProgBase.setquadobj!(m,diagm(0 => [10.0,10.0,10.0]))
            rows = [1, 2, 2, 2, 3, 3, 3]
            cols = [1, 1, 1, 2, 2, 3, 3]
            vals = Float64[2, 0.5, 0.5, 2, 1, 1, 1]
            MathProgBase.setquadobj!(m,rows,cols,vals)
            MathProgBase.optimize!(m)
            stat = MathProgBase.status(m)
            @test stat == :Optimal
            @test isapprox(MathProgBase.getobjval(m), 130/70, atol=1e-6)
            @test isapprox(norm(MathProgBase.getsolution(m) - [0.5714285714285715,0.4285714285714285,0.8571428571428572]), 0.0, atol=1e-6)
        end

        @testset "QP2" begin
            m = MathProgBase.LinearQuadraticModel(solver)
            MathProgBase.loadproblem!(m, [-1. 1.; 1. 1.], [0.,0.], [Inf,Inf], [1.,1.], [0.,0.], [Inf,Inf], :Max)
            MathProgBase.addquadconstr!(m, [2], [1.], [1], [1], [1.], '<', 2)
            MathProgBase.optimize!(m)
            stat = MathProgBase.status(m)
            @test stat == :Optimal
            @test isapprox(MathProgBase.getobjval(m), 2.25, atol=1e-6)
            @test isapprox(norm(MathProgBase.getsolution(m) - [0.5,1.75]), 0.0, atol=1e-3)
        end
    end
end

function qpdualtest(solver)
    @testset "Testing QP duals with $solver" begin
        # max x
        # s.t. x^2 <= 2
        @testset "QP1" begin
            m = MathProgBase.LinearQuadraticModel(solver)
            MathProgBase.loadproblem!(m, zeros(0,1), [-Inf], [Inf], [1.0], Float64[], Float64[], :Max)
            MathProgBase.addquadconstr!(m, [], [], [1], [1], [1.0], '<', 2.0)
            MathProgBase.optimize!(m)
            stat = MathProgBase.status(m)

            @test MathProgBase.numlinconstr(m) == 0
            @test MathProgBase.numquadconstr(m) == 1
            @test MathProgBase.numconstr(m) == 1
            @test stat == :Optimal
            @test isapprox(MathProgBase.getobjval(m), sqrt(2), atol=1e-6)
            @test isapprox(MathProgBase.getsolution(m)[1], sqrt(2), atol=1e-6)
            @test isapprox(MathProgBase.getquadconstrduals(m)[1], 0.5/sqrt(2), atol=1e-6)
        end

        # min -x
        # s.t. x^2 <= 2
        @testset "QP2" begin
            m = MathProgBase.LinearQuadraticModel(solver)
            MathProgBase.loadproblem!(m, zeros(0,1), [-Inf], [Inf], [-1.0], Float64[], Float64[], :Min)
            MathProgBase.addquadconstr!(m, [], [], [1], [1], [1.0], '<', 2.0)
            MathProgBase.optimize!(m)
            stat = MathProgBase.status(m)

            @test MathProgBase.numlinconstr(m) == 0
            @test MathProgBase.numquadconstr(m) == 1
            @test MathProgBase.numconstr(m) == 1
            @test stat == :Optimal
            @test isapprox(MathProgBase.getobjval(m), -sqrt(2), atol=1e-6)
            @test isapprox(MathProgBase.getsolution(m)[1], sqrt(2), atol=1e-6)
            @test isapprox(MathProgBase.getquadconstrduals(m)[1], -0.5/sqrt(2), atol=1e-6)
        end
    end
end

function socptest(solver)
    @testset "Testing SOCP interface with $solver" begin
        # min t
        # s.t. x + y >= 1
        #      x^2 + y^2 <= t^2
        #      t >= 0
        m = MathProgBase.LinearQuadraticModel(solver)
        MathProgBase.loadproblem!(m, [ 1.0 1.0 0.0 ], [-Inf,-Inf,0.0], [Inf,Inf,Inf], [0.0,0.0,1.0], [1.0],[Inf], :Min)
        MathProgBase.addquadconstr!(m, [], [], [1,2,3], [1,2,3], [1.0,1.0,-1.0],'<',0.0)
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)

        @test stat == :Optimal
        @test isapprox(MathProgBase.getobjval(m), sqrt(1/2), atol=1e-6)
        @test isapprox(norm(MathProgBase.getsolution(m) - [0.5,0.5,sqrt(1/2)]), 0.0, atol=1e-3)
    end
end
