using MathProgBase
using Compat
using Compat.Test

# Here the type represents the complete instance, but it
# could also store instance data.
mutable struct HS071 <: MathProgBase.AbstractNLPEvaluator
end

# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

function MathProgBase.initialize(d::HS071, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
end

MathProgBase.features_available(d::HS071) = [:Grad, :Jac, :Hess]

MathProgBase.eval_f(d::HS071, x) = x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]

function MathProgBase.eval_g(d::HS071, g, x)
    g[1] = x[1]   * x[2]   * x[3]   * x[4]
    g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

function MathProgBase.eval_grad_f(d::HS071, grad_f, x)
    grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
    grad_f[2] = x[1] * x[4]
    grad_f[3] = x[1] * x[4] + 1
    grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end

MathProgBase.jac_structure(d::HS071) = [1,1,1,1,2,2,2,2],[1,2,3,4,1,2,3,4]
# lower triangle only
MathProgBase.hesslag_structure(d::HS071) = [1,2,2,3,3,3,4,4,4,4],[1,1,2,1,2,3,1,2,3,4]


function MathProgBase.eval_jac_g(d::HS071, J, x)
    # Constraint (row) 1
    J[1] = x[2]*x[3]*x[4]  # 1,1
    J[2] = x[1]*x[3]*x[4]  # 1,2
    J[3] = x[1]*x[2]*x[4]  # 1,3
    J[4] = x[1]*x[2]*x[3]  # 1,4
    # Constraint (row) 2
    J[5] = 2*x[1]  # 2,1
    J[6] = 2*x[2]  # 2,2
    J[7] = 2*x[3]  # 2,3
    J[8] = 2*x[4]  # 2,4
end

function MathProgBase.eval_hesslag(d::HS071, H, x, σ, μ)
    # Again, only lower left triangle
    # Objective
    H[1] = σ * (2*x[4])               # 1,1
    H[2] = σ * (  x[4])               # 2,1
    H[3] = 0                          # 2,2
    H[4] = σ * (  x[4])               # 3,1
    H[5] = 0                          # 3,2
    H[6] = 0                          # 3,3
    H[7] = σ* (2*x[1] + x[2] + x[3])  # 4,1
    H[8] = σ * (  x[1])               # 4,2
    H[9] = σ * (  x[1])               # 4,3
    H[10] = 0                         # 4,4

    # First constraint
    H[2] += μ[1] * (x[3] * x[4])  # 2,1
    H[4] += μ[1] * (x[2] * x[4])  # 3,1
    H[5] += μ[1] * (x[1] * x[4])  # 3,2
    H[7] += μ[1] * (x[2] * x[3])  # 4,1
    H[8] += μ[1] * (x[1] * x[3])  # 4,2
    H[9] += μ[1] * (x[1] * x[2])  # 4,3

    # Second constraint
    H[1]  += μ[2] * 2  # 1,1
    H[3]  += μ[2] * 2  # 2,2
    H[6]  += μ[2] * 2  # 3,3
    H[10] += μ[2] * 2  # 4,4

end

function nlptest(solver)
    @testset "Testing NLP with $solver" begin
        m = MathProgBase.NonlinearModel(solver)
        l = [1,1,1,1]
        u = [5,5,5,5]
        lb = [25, 40]
        ub = [Inf, 40]
        MathProgBase.loadproblem!(m, 4, 2, l, u, lb, ub, :Min, HS071())
        MathProgBase.setwarmstart!(m,[1,5,5,1])

        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)

        @test stat == :Optimal
        x = MathProgBase.getsolution(m)
        @test isapprox(x[1], 1.0000000000000000, atol=1e-5)
        @test isapprox(x[2], 4.7429996418092970, atol=1e-5)
        @test isapprox(x[3], 3.8211499817883077, atol=1e-5)
        @test isapprox(x[4], 1.3794082897556983, atol=1e-5)
        @test isapprox(MathProgBase.getobjval(m), 17.014017145179164, atol=1e-5)

        # Test that a second call to optimize! works
        MathProgBase.setwarmstart!(m,[1,5,5,1])
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)
        @test stat == :Optimal
        if hasmethod(MathProgBase.freemodel!, Tuple{typeof(m)})
            MathProgBase.freemodel!(m)
        end
    end
end

# Same as above but no hessian callback
mutable struct HS071_2 <: MathProgBase.AbstractNLPEvaluator
end

# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

function MathProgBase.initialize(d::HS071_2, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
end

MathProgBase.features_available(d::HS071_2) = [:Grad, :Jac]

MathProgBase.eval_f(d::HS071_2, x) = MathProgBase.eval_f(HS071(), x)

MathProgBase.eval_g(d::HS071_2, g, x) = MathProgBase.eval_g(HS071(), g, x)

MathProgBase.eval_grad_f(d::HS071_2, grad_f, x) = MathProgBase.eval_grad_f(HS071(), grad_f, x)

MathProgBase.jac_structure(d::HS071_2) = MathProgBase.jac_structure(HS071())

MathProgBase.eval_jac_g(d::HS071_2, J, x) = MathProgBase.eval_jac_g(HS071(), J, x)

function nlptest_nohessian(solver)
    @testset "Testing NLP without Hessian with $solver" begin
        m = MathProgBase.NonlinearModel(solver)
        l = [1,1,1,1]
        u = [5,5,5,5]
        lb = [25, 40]
        ub = [Inf, 40]
        MathProgBase.loadproblem!(m, 4, 2, l, u, lb, ub, :Min, HS071_2())
        MathProgBase.setwarmstart!(m,[1,5,5,1])

        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)

        @test stat == :Optimal
        x = MathProgBase.getsolution(m)
        @test isapprox(x[1], 1.0000000000000000, atol=1e-5)
        @test isapprox(x[2], 4.7429996418092970, atol=1e-5)
        @test isapprox(x[3], 3.8211499817883077, atol=1e-5)
        @test isapprox(x[4], 1.3794082897556983, atol=1e-5)
        @test isapprox(MathProgBase.getobjval(m), 17.014017145179164, atol=1e-5)

        # Test that a second call to optimize! works
        MathProgBase.setwarmstart!(m,[1,5,5,1])
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)
        @test stat == :Optimal
        if hasmethod(MathProgBase.freemodel!, Tuple{typeof(m)})
            MathProgBase.freemodel!(m)
        end
    end
end

# a test for convex nonlinear solvers

mutable struct QCQP <: MathProgBase.AbstractNLPEvaluator
end

# min x - y
# st  x + x^2 + x*y + y^2 <= 1
#     -2 <= x, y <= 2
# solution: x+y = -1/3
# optimal objective -1-4/sqrt(3)

function MathProgBase.initialize(d::QCQP, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Jac, :Hess])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
end

MathProgBase.features_available(d::QCQP) = [:Grad, :Jac, :Hess]

MathProgBase.eval_f(d::QCQP, x) = x[1]-x[2]

function MathProgBase.eval_g(d::QCQP, g, x)
    g[1] = x[1] + x[1]^2 + x[1]*x[2] + x[2]^2
end

function MathProgBase.eval_grad_f(d::QCQP, grad_f, x)
    grad_f[1] = 1
    grad_f[2] = -1
end

MathProgBase.jac_structure(d::QCQP) = [1,1],[1,2]
# lower triangle only
MathProgBase.hesslag_structure(d::QCQP) = [1,2,2],[1,1,2]


function MathProgBase.eval_jac_g(d::QCQP, J, x)
    J[1] = 1 + 2x[1] + x[2]
    J[2] = x[1]+2x[2]
end

function MathProgBase.eval_hesslag(d::QCQP, H, x, σ, μ)
    # Again, only lower left triangle
    # Objective, linear

    # First constraint
    H[1] = 2μ[1] # d/dx^2
    H[2] = μ[1]  # d/dxdy
    H[3] = 2μ[1] # d/dy^2

end

function convexnlptest(solver)
    @testset "Testing convex NLP with $solver" begin
        m = MathProgBase.NonlinearModel(solver)
        l = [-2,-2]
        u = [2,2]
        lb = [-Inf]
        ub = [1.0]
        MathProgBase.loadproblem!(m, 2, 1, l, u, lb, ub, :Min, QCQP())

        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)

        @test stat == :Optimal
        x = MathProgBase.getsolution(m)
        @test isapprox(x[1]+x[2], -1/3, atol=1e-3)
        @test isapprox(MathProgBase.getobjval(m), -1-4/sqrt(3), atol=1e-5)

        # Test that a second call to optimize! works
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)
        @test stat == :Optimal
        if hasmethod(MathProgBase.freemodel!, Tuple{typeof(m)})
            MathProgBase.freemodel!(m)
        end
    end
end


# a test for unconstrained nonlinear solvers

mutable struct Rosenbrock <: MathProgBase.AbstractNLPEvaluator
end

# min (1-x)^2 + 100*(y-x^2)^2
# solution: (x,y) = (1,1)
# optimal objective 0

function MathProgBase.initialize(d::Rosenbrock, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad, :Hess])
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
end

MathProgBase.features_available(d::Rosenbrock) = [:Grad, :Hess]

MathProgBase.eval_f(d::Rosenbrock, x) = (1-x[1])^2 + 100*(x[2]-x[1]^2)^2

MathProgBase.eval_g(d::Rosenbrock, g, x) = nothing

function MathProgBase.eval_grad_f(d::Rosenbrock, grad_f, x)
    grad_f[1] = -2*(1-x[1]) - 400*x[1]*(x[2]-x[1]^2)
    grad_f[2] = 200*(x[2]-x[1]^2)
end

MathProgBase.hesslag_structure(d::Rosenbrock) = [1,2,2],[1,1,2]

function MathProgBase.eval_hesslag(d::Rosenbrock, H, x, σ, μ)
    # Again, only lower left triangle

    H[1] = σ*(2+100*(8x[1]^2-4*(x[2]-x[1]^2))) # d/dx^2
    H[2] = -σ*400*x[1]  # d/dxdy
    H[3] = σ*200 # d/dy^2

end

function rosenbrocktest(solver)
    @testset "Testing NLP on the Rosenbrock function with $solver" begin
        m = MathProgBase.NonlinearModel(solver)
        l = [-Inf,-Inf]
        u = [Inf,Inf]
        lb = Float64[]
        ub = Float64[]
        MathProgBase.loadproblem!(m, 2, 0, l, u, lb, ub, :Min, Rosenbrock())

        MathProgBase.setwarmstart!(m,[10.0,10.0])
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)

        @test stat == :Optimal
        x = MathProgBase.getsolution(m)
        @test isapprox(x[1], 1.0, atol=1e-5)
        @test isapprox(x[2], 1.0, atol=1e-5)
        @test isapprox(MathProgBase.getobjval(m), 0.0, atol=1e-5)

        # Test that a second call to optimize! works
        MathProgBase.optimize!(m)
        stat = MathProgBase.status(m)
        @test stat == :Optimal
        if hasmethod(MathProgBase.freemodel!, Tuple{typeof(m)})
            MathProgBase.freemodel!(m)
        end
    end
end
