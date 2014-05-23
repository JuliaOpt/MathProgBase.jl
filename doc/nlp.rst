---------------------
Nonlinear Programming
---------------------

MathProgBase provides an interface for nonlinear programming which is independent of both the solver and the user's representation of the problem, whether using an algebraic modeling language or customized low-level code.

The diagram below illustrates MathProgBase as the connection between typical NLP solvers IPOPT, MOSEK, and KNITRO, and modeling languages such as `JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ and `AMPL <http://ampl.com/>`_ (via `ampl.jl <https://github.com/dpo/ampl.jl>`_).

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

Where :math:`x \in \mathbb{R}^n, f: \mathbb{R}^n \to \mathbb{R}, g: \mathbb{R}^n \to \mathbb{R}^m`, and vectors :math:`lb \in \mathbb{R}^m \cup \{-\infty\}, ub \in \mathbb{R}^m \cup \{\infty\},l \in \mathbb{R}^n \cup \{-\infty\}, u \in \mathbb{R}^n \cup \{\infty\}`.

The objective function :math:`f` and constraint function :math:`g` may be nonlinear and nonconvex, but are typically expected to be twice differentiable.

Below we describe extensions to the ``MathProgSolverInterface`` for these nonlinear programming problems.

.. function:: loadnonlinearproblem!(m::AbstractMathProgModel, numVar, numConstr, l, u, lb, ub, d::AbstractNLPEvaluator)
    
    Loads the nonlinear programming problem into the model. The parameter `numVar` is the number of variables in the problem, ``numConstr`` is the number of constraints, ``l`` contains the variable lower bounds, ``u`` contains the variable upper bounds, ``lb`` contains the constraint lower bounds, and ``ub`` contains the constraint upper bounds. The final parameter ``d`` is an instance of an ``AbstractNLPEvaluator``, described below, which may be queried for evaluating :math:`f` and :math:`g` and their corresponding derivatives.

The abstract type ``AbstractNLPEvaluator`` is used by solvers for accessing the objective function :math:`f` and constraints :math:`g`. Solvers may query the value, gradients, Hessian-vector products, and the Hessian of the Lagrangian.

.. function:: initialize(d::AbstractNLPEvaluator, requested_features::Vector{Symbol})

    Must be called before any other methods. The vector ``requested_features``
    lists features requested by the solver. These may include ``:Grad`` for gradients
    of :math:`f` and Jacobians of :math:`g`, ``:HessVec`` for Hessian-vector
    and Hessian-of-Lagrangian-vector products, ``:Hess`` for full Hessians and
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
    structurally nonzero element. These indices may not be sorted and can contain
    duplicates, in which case the solver should combine the corresponding elements by
    adding them together.

.. function:: hesslag_structure(d::AbstractNLPEvaluator)

    Returns the sparsity structure of the Hessian-of-the-Lagrangian matrix 
    :math:`\nabla^2 f + \sum_{i=1}^m \nabla^2 g_i` as a tuple ``(I,J)``
    where ``I`` contains the row indices and ``J`` contains the column indices of each
    structurally nonzero element. These indices may not be sorted and can contain
    duplicates, in which case the solver should combine the corresponding elements by
    adding them together. Any mix of lower and upper-triangular indices is valid.
    Elements ``(i,j)`` and ``(j,i)``, if both present, should be treated as duplicates.

.. function:: eval_jac_g(d::AbstractNLPEvaluator, J, x)

    Evaluates the sparse Jacobian matrix :math:`J_g(x) = \left[ \begin{array}{c} \nabla g_1(x) \\ \nabla g_2(x) \\ \vdots \\ \nabla g_m(x) \end{array}\right]`.
    The result is stored in the vector ``J`` in the same order as the indices returned
    by ``jac_structure``.

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

.. function:: obj_expr(d::AbstractNLPEvaluator)

    Returns an expression graph for the objective function. *FORMAT TO BE DETERMINED*

.. function:: constr_expr(d::AbstractNLPEvaluator, i)

    Returns an expression graph for the :math:`i\text{th}` constraint. *FORMAT TO BE DETERMINED*


The solution vector, optimal objective value, termination status, etc. should be accessible from the standard methods, e.g., ``getsolution``, ``getobjval``, ``status``, respectively.

