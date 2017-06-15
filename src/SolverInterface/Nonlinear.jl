# Methods for the Nonlinear interface

abstract type AbstractNonlinearModel <: AbstractMathProgModel end
export AbstractNonlinearModel

abstract type AbstractNLPEvaluator end
export AbstractNLPEvaluator


###  NLP Attributes ###

"""
    ConstraintNLPDual(N)
    ConstraintNLPDual()

The assignment to the NLP constraint dual values in result `N`. If `N` is omitted, it is 1 by default.
"""
struct ConstraintNLPDual <: AbstractAttribute
    N::Int
end
ConstraintNLPDual() = ConstraintNLPDual(1)


"""
    ConstraintNLPDualStart()

An initial assignment of the NLP constriant duals that the solver may use
to warm-start the solve.
"""
struct ConstraintNLPDualStart <: AbstractAttribute end


### methods for AbstractNLPEvaluator ###

"""

"""
function NonlinearModel end


"""
    loadnlp!(m::AbstractNonlinearModel, numVar, numConstr, l, u, lb, ub, sense, d::AbstractNLPEvaluator)

Loads the nonlinear programming problem into the model. The parameter `numVar` is the number of variables in the problem, `numConstr` is the number of constraints, `l` contains the variable lower bounds, `u` contains the variable upper bounds, `lb` contains the constraint lower bounds, and `ub` contains the constraint upper bounds. Sense contains the symbol `:Max` or `:Min`, indicating the direction of optimization. The final parameter `d` is an instance of an `AbstractNLPEvaluator`, described below, which may be queried for evaluating ``f`` and ``g`` and their corresponding derivatives.
"""
function loadnlp! end


"""
    initialize(d::AbstractNLPEvaluator, requested_features::Vector{Symbol})

Must be called before any other methods. The vector `requested_features`
lists features requested by the solver. These may include `:Grad` for gradients
of ``f``, `:Jac` for explicit Jacobians of ``g``, `:JacVec` for
Jacobian-vector products, `:HessVe` for Hessian-vector
and Hessian-of-Lagrangian-vector products, `:Hess` for explicit Hessians and
Hessian-of-Lagrangians, and `:ExprGraph` for expression graphs.
"""
function initialize end


"""
    features_available(d::AbstractNLPEvaluator)

Returns the subset of features available for this problem instance, as a
list of symbols in the same format as in `initialize`.
"""
function features_available end


"""
    eval_f(d::AbstractNLPEvaluator, x)

Evaluate ``f(x)``, returning a scalar value.
"""
function eval_f end


"""
    eval_g(d::AbstractNLPEvaluator, g, x)

Evaluate ``g(x)``, storing the result in the vector `g` which
must be of the appropriate size.
"""
function eval_g end


"""
    eval_grad_f(d::AbstractNLPEvaluator, g, x)

Evaluate ``\\nabla f(x)`` as a dense vector, storing 
the result in the vector `g` which must be of the appropriate size.
"""
function eval_grad_f end


"""
    jac_structure(d::AbstractNLPEvaluator)

Returns the sparsity structure of the Jacobian matrix, ``J_g(x) = \\left[ \\begin{array}{c} \\nabla g_1(x) \\\\ \\nabla g_2(x) \\\\ \\vdots \\\\ \\nabla g_m(x) \\end{array} \\right]`` where ``g_i`` is the ``i\\text{th}`` component of ``g``. The sparsity structure is assumed to be independent of the point ``x``. Returns a tuple ``(I,J)`` where ``I`` contains the row indices and ``J`` contains the column indices of each structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together.
"""
function jac_structure end


"""
    hesslag_structure(d::AbstractNLPEvaluator)

Returns the sparsity structure of the Hessian-of-the-Lagrangian matrix 
``\\nabla^2 f + \\sum_{i=1}^m \\nabla^2 g_i`` as a tuple ``(I,J)``
where ``I`` contains the row indices and ``J`` contains the column indices of each
structurally nonzero element. These indices are not required to be sorted and can contain duplicates, in which case the solver should combine the corresponding elements by adding them together. Any mix of lower and upper-triangular indices is valid.
Elements ``(i,j)`` and ``(j,i)``, if both present, should be treated as duplicates.
"""
function hesslag_structure end


"""
    eval_jac_g(d::AbstractNLPEvaluator, J, x)

Evaluates the sparse Jacobian matrix ``J_g(x) = \\left[ \\begin{array}{c} \\nabla g_1(x) \\\\ \\nabla g_2(x) \\\\ \\vdots \\\\ \\nabla g_m(x) \\end{array} \\right]``.
The result is stored in the vector `J` in the same order as the indices returned
by `jac_structure`.
"""
function eval_jac_g end


"""
    eval_jac_prod(d::AbstractNLPEvaluator, y, x, w)

Computes the Jacobian-vector product ``J_g(x)w``, storing the result in the vector `y`.
"""
function eval_jac_prod end


"""
    eval_jac_prod_t(d::AbstractNLPEvaluator, y, x, w)

Computes the Jacobian-transpose-vector product ``J_g(x)^T w``,
storing the result in the vector `y`.

"""
function eval_jac_prod_t end


"""
    eval_hesslag_prod(d::AbstractNLPEvaluator, h, x, v, σ, μ)

Given scalar weight ``σ`` and vector of constraint weights ``μ``, computes the Hessian-of-the-Lagrangian-vector product ``\\left( \\sigma \\nabla^2 f(x) + \\sum_{i=1}^m \\mu_i \\nabla^2 g_i(x) \\right)v``, 
storing the result in the vector ``h``.
"""
function eval_hesslag_prod end


"""
    eval_hesslag(d::AbstractNLPEvaluator, H, x, σ, μ)

Given scalar weight `σ` and vector of constraint weights `μ`, 
computes the sparse Hessian-of-the-Lagrangian matrix 
``\\sigma \\nabla^2 f(x) + \\sum_{i=1}^m \\mu_i \\nabla^2 g_i(x)``, 
storing the result in the vector `H` in the same order as the indices
returned by `hesslag_structure`.
"""
function eval_hesslag end


"""
    obj_expr(d::AbstractNLPEvaluator)

Returns an expression graph for the objective function as a standard Julia `Expr`
object. All sums and products are flattened out as simple `Expr(:+,...)` and
`Expr(:*,...)` objects. The symbol `x` is used as a placeholder for the
vector of decision variables. No other undefined symbols are permitted;
coefficients are embedded as explicit values.
For example, the expression
``x_1+\\sin(x_2/\\exp(x_3))`` would be represented as the Julia object
`:(x[1] + sin(x[2]/exp(x[3])))`. See the [Julia manual](http://docs.julialang.org/en/release-0.3/manual/metaprogramming/#expressions-and-eval) for more information
on the structure of `Expr` objects. There are currently no restrictions on
recognized functions; typically these will be built-in Julia functions like
`^`, `exp`, `log`, `cos`, `tan`, `sqrt`, etc., but modeling
interfaces may choose to extend these basic functions.
"""
function obj_expr end


"""
    constr_expr(d::AbstractNLPEvaluator, i)

Returns an expression graph for the ``i\\text{th}`` constraint in the same format as described above. The head of the expression is ``:comparison``, indicating the sense
of the constraint. The right-hand side of the comparison must be a constant; that is,
`:(x[1]^3 <= 1)` is allowed, while `:(1 <= x[1]^3)` is not valid.
Double-sided constraints are allowed, in which case both the lower bound and
upper bounds should be constants; for example, `:(-1 <= cos(x[1]) + sin(x[2]) <= 1)` is valid.
"""
function constr_expr end


### fallback options ###
"""
    isobjlinear(::AbstractNLPEvaluator)

`true` if the objective function is known to be linear,
`false` otherwise.
"""
isobjlinear(::AbstractNLPEvaluator) = false


"""
    isobjquadratic(::AbstractNLPEvaluator)

`true` if the objective function is known to be quadratic (convex or nonconvex), `false` otherwise.
"""
isobjquadratic(::AbstractNLPEvaluator) = false


"""
    isconstrlinear(::AbstractNLPEvaluator, i::Integer)

`true` if the ``i\text{th}`` constraint is known to be linear, `false` otherwise.
"""
isconstrlinear(::AbstractNLPEvaluator, i::Integer) = false



