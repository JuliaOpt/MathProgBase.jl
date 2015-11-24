----------------
Nonlinear models
----------------

The abstraction for nonlinear models is independent of both the solver and the user's representation of the problem, whether using an algebraic modeling language or customized low-level code.

The diagram below illustrates MathProgBase as the connection between typical nonlinear programming (NLP) solvers IPOPT, MOSEK, and KNITRO, and modeling languages such as `JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ and `AMPL <http://ampl.com/>`_ (via `AmplNLReader.jl <https://github.com/dpo/AmplNLReader.jl>`_).

.. graph:: foo
  
   node [shape="box"];

   subgraph clusterA {
   
    "IPOPT" -- "MOSEK" -- "KNITRO" -- "..." [style="invis", constraint="false"];
    label="Solvers";
    penwidth=0;

    }
    "IPOPT" -- "JuMP" [style="invis"];
   "IPOPT" -- "MathProgBase";
   "MOSEK" -- "MathProgBase";
   "KNITRO" -- "MathProgBase";
   "..." -- "MathProgBase";
   "MathProgBase" -- "JuMP";
   "MathProgBase" -- "AMPL";
   "MathProgBase" -- "User";
   "MathProgBase" -- x;

    subgraph clusterB {
    x [label="..."];
    "JuMP" -- "AMPL" -- "User" -- x [style="invis", constraint="false"];
    label="Modeling";
    penwidth=0;
    }
    rankdir=LR;
 
This structure also makes it easy to connect solvers written in Julia itself with user-provided instances in a variety of formats.

We take the prototypical format for a nonlinear problem as

.. math::
    \min_{x}\, &f(x)\\
    s.t.     &lb \leq g(x) \leq ub\\
             &l \leq x \leq u\\

Where :math:`x \in \mathbb{R}^n, f: \mathbb{R}^n \to \mathbb{R}, g: \mathbb{R}^n \to \mathbb{R}^m`, and vectors :math:`lb \in (\mathbb{R} \cup \{-\infty\})^m, ub \in (\mathbb{R} \cup \{\infty\})^m,l \in (\mathbb{R} \cup \{-\infty\})^n, u \in (\mathbb{R} \cup \{\infty\})^n`.

The objective function :math:`f` and constraint function :math:`g` may be nonlinear and nonconvex, but are typically expected to be twice differentiable.

Below we describe the interface for ``AbstractNonlinearModel``.

.. function:: loadproblem!(m::AbstractNonlinearModel, numVar, numConstr, l, u, lb, ub, sense, d::AbstractNLPEvaluator)
    
    Loads the nonlinear programming problem into the model. The parameter `numVar` is the number of variables in the problem, ``numConstr`` is the number of constraints, ``l`` contains the variable lower bounds, ``u`` contains the variable upper bounds, ``lb`` contains the constraint lower bounds, and ``ub`` contains the constraint upper bounds. Sense contains the symbol ``:Max`` or ``:Min``, indicating the direction of optimization. The final parameter ``d`` is an instance of an ``AbstractNLPEvaluator``, described below, which may be queried for evaluating :math:`f` and :math:`g` and their corresponding derivatives.

The abstract type ``AbstractNLPEvaluator`` is used by solvers for accessing the objective function :math:`f` and constraints :math:`g`. Solvers may query the value, gradients, Hessian-vector products, and the Hessian of the Lagrangian.

.. function:: initialize(d::AbstractNLPEvaluator, requested_features::Vector{Symbol})

    Must be called before any other methods. The vector ``requested_features``
    lists features requested by the solver. These may include ``:Grad`` for gradients
    of :math:`f`, ``:Jac`` for explicit Jacobians of :math:`g`, ``:JacVec`` for
    Jacobian-vector products, ``:HessVec`` for Hessian-vector
    and Hessian-of-Lagrangian-vector products, ``:Hess`` for explicit Hessians and
    Hessian-of-Lagrangians, and ``:ExprGraph`` for expression graphs.

.. function:: features_available(d::AbstractNLPEvaluator)

    Returns the subset of features available for this problem instance, as a
    list of symbols in the same format as in ``initialize``.

.. function:: eval_f(d::AbstractNLPEvaluator, x)

    Evaluate :math:`f(x)`, returning a scalar value.

.. function:: eval_g(d::AbstractNLPEvaluator, g, x)

    Evaluate :math:`g(x)`, storing the result in the vector ``g`` which
    must be of the appropriate size.

.. function:: eval_grad_f(d::AbstractNLPEvaluator, g, x)

    Evaluate :math:`\nabla f(x)` as a dense vector, storing 
    the result in the vector ``g`` which must be of the appropriate size.

.. function:: jac_structure(d::AbstractNLPEvaluator)

    Returns the sparsity structure of the Jacobian matrix :math:`J_g(x) = \left[ \begin{array}{c} \nabla g_1(x) \\ \nabla g_2(x) \\ \vdots \\ \nabla g_m(x) \end{array}\right]` where :math:`g_i` is the :math:`i\text{th}` component of :math:`g`. The sparsity structure
    is assumed to be independent of the point :math:`x`. Returns a tuple ``(I,J)``
    where ``I`` contains the row indices and ``J`` contains the column indices of each
    structurally nonzero element. These indices are not required to be sorted and can contain
    duplicates, in which case the solver should combine the corresponding elements by
    adding them together.

.. function:: hesslag_structure(d::AbstractNLPEvaluator)

    Returns the sparsity structure of the Hessian-of-the-Lagrangian matrix 
    :math:`\nabla^2 f + \sum_{i=1}^m \nabla^2 g_i` as a tuple ``(I,J)``
    where ``I`` contains the row indices and ``J`` contains the column indices of each
    structurally nonzero element. These indices are not required to be sorted and can contain
    duplicates, in which case the solver should combine the corresponding elements by
    adding them together. Any mix of lower and upper-triangular indices is valid.
    Elements ``(i,j)`` and ``(j,i)``, if both present, should be treated as duplicates.

.. function:: eval_jac_g(d::AbstractNLPEvaluator, J, x)

    Evaluates the sparse Jacobian matrix :math:`J_g(x) = \left[ \begin{array}{c} \nabla g_1(x) \\ \nabla g_2(x) \\ \vdots \\ \nabla g_m(x) \end{array}\right]`.
    The result is stored in the vector ``J`` in the same order as the indices returned
    by ``jac_structure``.

.. function:: eval_jac_prod(d::AbstractNLPEvaluator, y, x, w)

    Computes the Jacobian-vector product :math:`J_g(x)w`,
    storing the result in the vector ``y``.

.. function:: eval_jac_prod_t(d::AbstractNLPEvaluator, y, x, w)

    Computes the Jacobian-transpose-vector product :math:`J_g(x)^Tw`,
    storing the result in the vector ``y``.

.. function:: eval_hesslag_prod(d::AbstractNLPEvaluator, h, x, v, σ, μ)

    Given scalar weight ``σ`` and vector of constraint weights ``μ``, 
    computes the Hessian-of-the-Lagrangian-vector product 
    :math:`\left(\sigma\nabla^2 f(x) + \sum_{i=1}^m \mu_i \nabla^2 g_i(x)\right)v`, 
    storing the result in the vector ``h``.

.. function:: eval_hesslag(d::AbstractNLPEvaluator, H, x, σ, μ)

    Given scalar weight ``σ`` and vector of constraint weights ``μ``, 
    computes the sparse Hessian-of-the-Lagrangian matrix 
    :math:`\sigma\nabla^2 f(x) + \sum_{i=1}^m \mu_i \nabla^2 g_i(x)`, 
    storing the result in the vector ``H`` in the same order as the indices
    returned by ``hesslag_structure``.

.. function:: isobjlinear(d::AbstractNLPEvaluator)

    ``true`` if the objective function is known to be linear,
    ``false`` otherwise.

.. function:: isobjquadratic(d::AbstractNLPEvaluator)

    ``true`` if the objective function is known to be quadratic (convex or nonconvex),
    ``false`` otherwise.

.. function:: isconstrlinear(d::AbstractNLPEvaluator, i)

    ``true`` if the :math:`i\text{th}` constraint is known to be linear,
    ``false`` otherwise.

.. function:: obj_expr(d::AbstractNLPEvaluator)

    Returns an expression graph for the objective function as a standard Julia ``Expr``
    object. All sums and products are flattened out as simple ``Expr(:+,...)`` and
    ``Expr(:*,...)`` objects. The symbol ``x`` is used as a placeholder for the
    vector of decision variables. No other undefined symbols are permitted;
    coefficients are embedded as explicit values.
    For example, the expression
    :math:`x_1+\sin(x_2/\exp(x_3))` would be represented as the Julia object
    ``:(x[1] + sin(x[2]/exp(x[3])))``. See the `Julia manual <http://docs.julialang.org/en/release-0.3/manual/metaprogramming/#expressions-and-eval>`_ for more information
    on the structure of ``Expr`` objects. There are currently no restrictions on
    recognized functions; typically these will be built-in Julia functions like
    ``^``, ``exp``, ``log``, ``cos``, ``tan``, ``sqrt``, etc., but modeling
    interfaces may choose to extend these basic functions.

.. function:: constr_expr(d::AbstractNLPEvaluator, i)

    Returns an expression graph for the :math:`i\text{th}` constraint in the same format as described above. The head of the expression is ``:comparison``, indicating the sense
    of the constraint. The right-hand side of the comparison must be a constant; that is,
    ``:(x[1]^3 <= 1)`` is allowed, while ``:(1 <= x[1]^3)`` is not valid.
    Double-sided constraints are allowed, in which case both the lower bound and
    upper bounds should be constants; for example, ``:(-1 <= cos(x[1]) + sin(x[2]) <= 1)`` is valid.

Nonlinear solvers may also provide optimal Lagrange multipliers if available through ``getreducedcosts`` and ``getconstrduals``.

.. function:: getreducedcosts(m::AbstractNonlinearModel)

    Returns the dual solution vector corresponding to the variable bounds,
    known as the reduced costs. Not available when integer variables are present.

.. function:: getconstrduals(m::AbstractNonlinearModel)

    Returns the dual solution vector corresponding to the constraints.
    Not available when integer variables are present.



