------------------
Linear Programming
------------------

.. function:: linprog(c, A, sense, b, l, u, solver)

Solves the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &a_i^Tx \text{ sense}_i \, b_i \forall\,\, i\\
             &l \leq x \leq u\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix, with rows :math:`a_i` (viewed as column-oriented vectors)
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` specifies the desired solver, see :ref:`choosing solvers <choosing-solvers>`.

A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. note::
    Linear programming solvers extensively exploit the sparsity of the constraint matrix ``A``. While both dense and sparse matrices are accepted, for large-scale problems sparse matrices should be provided if permitted by the problem structure.

A shortened version is defined as::

    linprog(c, A, sense, b, solver) = linprog(c, A, sense, b, 0, Inf, solver)

The ``linprog`` function returns an instance of the type::

    type LinprogSolution
        status
        objval
        sol
        attrs
    end

where ``status`` is a termination status symbol, one of ``:Optimal``, ``:Infeasible``, ``:Unbounded``, ``:UserLimit`` (iteration limit or timeout), ``:Error`` (and maybe others).

If ``status`` is ``:Optimal``, the other members have the following values:

* ``objval`` -- optimal objective value
* ``sol`` -- primal solution vector
* ``attrs`` -- a dictionary that may contain other relevant attributes such as:

  - ``redcost`` -- dual multipliers for active variable bounds (zero if inactive)
  - ``lambda`` -- dual multipliers for active linear constraints (equalities are always active)

If ``status`` is ``:Infeasible``, the ``attrs`` member will contain an ``infeasibilityray`` if available; similarly for ``:Unbounded`` problems, ``attrs`` will contain an ``unboundedray`` if available.

..
  - ``colbasis`` -- optimal simplex basis statuses for the variables (columns) if available. Possible values are ``:NonbasicAtLower``, ``:NonbasicAtUpper``, ``:Basic``, and ``:Superbasic`` (not yet implemented by any solvers)
  - ``rowbasis`` -- optimal simplex basis statuses for the constraints (rows) if available (not yet implemented by any solvers)

For example, we can solve the two-dimensional problem (see ``test/linprog.jl``):

.. math::
    \min_{x,y}\, &-x\\
    s.t.         &2x + y \leq 1.5\\
                 & x \geq 0, y \geq 0

by::

    using MathProgBase, Clp

    sol = linprog([-1,0],[2 1],'<',1.5, ClpSolver())
    if sol.status == :Optimal
        println("Optimal objective value is $(sol.objval)")
        println("Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2])]")
    else
        println("Error: solution status $(sol.status)")
    end

.. function:: linprog(c, A, lb, ub, l, u, solver)

This variant allows one to specify two-sided linear constraints (also known as range constraints)
to solve the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &lb \leq Ax \leq ub\\
             &l \leq x \leq u\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix
*    ``lb`` is the vector of row lower bounds
*    ``ub`` is the vector of row upper bounds
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` specifies the desired solver, see :ref:`choosing solvers <choosing-solvers>`.

A scalar is accepted for the ``l``, ``u``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

A shortened version is defined as::

    linprog(c, A, lb, ub, solver) = linprog(c, A, lb, ub, 0, Inf, solver)

.. note::
    The function ``linprog`` calls two independent functions for building and solving the linear programming problem, namely ``buildlp`` and ``solvelp``.

.. function:: buildlp(c, A, sense, b, l, u, solver)

Builds the linear programming problem as defined in ``linprog`` and accepts the following arguments:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` specifies the desired solver, see :ref:`choosing solvers <choosing-solvers>`.

A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. function:: buildlp(c, A, lb, ub, l, u, solver)

This variant of ``buildlp`` allows to specify two-sided linear constraints (also known as range constraints) similar to ``linprog``, and accepts the following arguments:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix
*    ``lb`` is the vector of row lower bounds
*    ``ub`` is the vector of row upper bounds
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` specifies the desired solver, see :ref:`choosing solvers <choosing-solvers>`.

A scalar is accepted for the ``l``, ``u``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

The ``buildlp`` function returns an ``AbstractLinearQuadraticModel`` that can be input to ```solvelp``` in order to obtain a solution.

.. function:: solvelp(m)

Solves the linear programming problem as defined in ``linprog``` and accepts the following argument:

*    ``m`` is an ``AbstractLinearQuadraticModel`` (e.g., as returned by ``buildlp``).

The ``solvelp`` function returns an instance of the type::

    type LinprogSolution
        status
        objval
        sol
        attrs
    end
