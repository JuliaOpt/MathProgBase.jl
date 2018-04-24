using Compat.Test
using Compat.LinearAlgebra
using Compat.SparseArrays
using MathProgBase

function linprogtest(solver; objtol = 1e-7, primaltol = 1e-6)
    @testset "Testing linprog and subfunctions with $solver" begin
        # min -x
        # s.t. 2x + y <= 1.5
        # x,y >= 0
        # solution is (0.75,0) with objval -0.75

        @testset "LP1" begin
            @testset "with buildlp and solvelp" begin
                m = buildlp([-1,0],[2 1],'<',1.5,solver)
                sol = solvelp(m)
                @test sol.status == :Optimal
                @test isapprox(sol.objval, -0.75, atol=objtol)
                @test isapprox(norm(sol.sol - [0.75,0.0]), 0.0, atol=primaltol)
            end

            @testset "with linprog" begin
                sol = linprog([-1,0],[2 1],'<',1.5,solver)
                @test sol.status == :Optimal
                @test isapprox(sol.objval, -0.75, atol=objtol)
                @test isapprox(norm(sol.sol - [0.75,0.0]), 0.0, atol=primaltol)
            end
        end

        @testset "LP2" begin
            @testset "with buildlp and solvelp" begin
                m = buildlp([-1,0],[2 1],'<',1.5,solver)
                sol = solvelp(m)
                @test sol.status == :Optimal
                @test isapprox(sol.objval, -0.75, atol=objtol)
                @test isapprox(norm(sol.sol - [0.75,0.0]), 0.0, atol=primaltol)
            end

            @testset "with linprog" begin
                sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
                @test sol.status == :Optimal
                @test isapprox(sol.objval, -0.75, atol=objtol)
                @test isapprox(norm(sol.sol - [0.75,0.0]), 0.0, atol=primaltol)
            end
        end

        # test infeasible problem:
        # min x
        # s.t. 2x+y <= -1
        # x,y >= 0

        @testset "LP3 infeasible" begin
            @testset "with buildlp and solvelp" begin
                m = buildlp([1,0],[2 1],'<',-1,solver)
                sol = solvelp(m)
                @test sol.status == :Infeasible

                r = sol.attrs[:infeasibilityray]
                @test length(r) == 1
                r[1] = sol.attrs[:infeasibilityray][1]
                @test isapprox(r[1]/abs(r[1]), -1.0)
            end

            @testset "with linprog" begin
                sol = linprog([1,0],[2 1],'<',-1,solver)
                @test sol.status == :Infeasible

                r = sol.attrs[:infeasibilityray]
                @test length(r) == 1
                r[1] = sol.attrs[:infeasibilityray][1]
                @test isapprox(r[1]/abs(r[1]), -1.0)
            end
        end

        # test unbounded problem:
        # min -x-y
        # s.t. -x+2y <= 0
        # x,y >= 0

        @testset "LP4 unbounded" begin
            @testset "with buildlp and solvelp" begin
                m = buildlp([-1,-1],[-1 2],'<',[0],solver)
                sol = solvelp(m)
                @test sol.status == :Unbounded
            end

            @testset "with linprog" begin
                sol = linprog([-1,-1],[-1 2],'<',[0],solver)
                @test sol.status == :Unbounded
            end
        end

        # unbounded problem with unique ray:
        # min -x-y
        # s.t. x-y == 0
        # x,y >= 0

        @testset "LP4 unbounded with unique ray" begin
            @testset "with buildlp and solvelp" begin
                m = buildlp([-1,-1],[1 -1],'=',0,solver)
                sol = solvelp(m)
                @test sol.status == :Unbounded
                @test sol.attrs[:unboundedray][1] > 1e-7
                @test isapprox(sol.attrs[:unboundedray][1], sol.attrs[:unboundedray][2])
            end

            @testset "with linprog" begin
                sol = linprog([-1,-1],[1 -1],'=',0,solver)
                @test sol.status == :Unbounded
                @test sol.attrs[:unboundedray][1] > 1e-7
                @test isapprox(sol.attrs[:unboundedray][1], sol.attrs[:unboundedray][2])
            end
        end
    end
end
