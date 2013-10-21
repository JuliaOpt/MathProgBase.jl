===============
MathProgBase.jl
===============

.. module:: MathProgBase
   :synopsis: Solver-independent functions and low-level interfaces for Mathematical Programming

`MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl>`_ provides high-level
one-shot functions for linear and mixed-integer programming,
as well as a solver-independent low-level interface for implementing advanced techniques
that require efficiently solving a sequence of linear programming problems.

To use MathProgBase, an external solver must be installed. See CHOOSING SOLVERS.

------------------
Linear Programming
------------------

.. function:: linprog(c, A, sense, b, lb, ub, solver)

Solves the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &a_i^Tx \text{ sense}_i \, b_i \forall\,\, i\\
             &lb \leq x \leq ub\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix, with rows :math:`a_i`
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``lb`` is the vector of lower bounds on the variables
*    ``ub`` is the vector of upper bounds on the variables, and
*    ``solver`` is an *optional* parameter specifying the desired solver, see CHOOSING SOLVERS. If this parameter is not provided, the default solver is used.
 
A scalar is accepted for the ``b``, ``sense``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. note::
    Linear programming solvers extensively exploit the sparsity of the constraint matrix ``A``. While both dense and sparse matrices are acceped, for large-scale problems sparse matrices should be provided if permitted by the problem structure.  

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
  - ``colbasis`` -- optimal simplex basis statuses for the variables (columns) if available. Possible values are ``:NonbasicAtLower``, ``:NonbasicAtUpper``, ``:Basic``, and ``:Superbasic`` (not yet implemented by any solvers)
  - ``rowbasis`` -- optimal simplex basis statuses for the constraints (rows) if available (not yet implemented by any solvers)

For example, we can solve the two-dimensional problem (see ``test/linprog.jl``):

.. math::
    \min_{x,y}\, &-x\\
    s.t.         &2x + y \leq 1.5\\
                 & x \geq 0, y \geq 0

by::

    using MathProgBase
    
    sol = linprog([-1,0],[2 1],'<',1.5)
    if sol.status == :Optimal
        println("Optimal objective value is $(sol.objval)")
        println("Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2])]")
    else
        println("Error: solution status $(sol.status)")
    end

.. function:: linprog(c, A, l, u, lb, ub, solver)

This variant allows one to specify two-sided linear constraints (also known as range constraints)
to solve the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &l \leq Ax \leq u\\
             &lb \leq x \leq ub\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``A`` is the constraint matrix
*    ``l`` is the vector of row lower bounds
*    ``u`` is the vector of row upper bounds
*    ``lb`` is the vector of lower bounds on the variables
*    ``ub`` is the vector of upper bounds on the variables, and
*    ``solver`` is an *optional* parameter specifying the desired solver, see CHOOSING SOLVERS. If this parameter is not provided, the default solver is used.
 
A scalar is accepted for the ``l``, ``u``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

A shortened version is defined as::

    linprog(c, A, l, u, solver) = linprog(c, A, l, u, 0, Inf, solver)


-------------------------
Mixed-integer Programming
-------------------------

.. function:: mixintprog(c, A, sense, b, vartypes, lb, ub, solver)

Solves the same optimization problem as ``linprog`` above, except variables
are additionally constrained to take only integer values if the corresponding
entry in the ``varypes`` vector is the character ``'I'``. Continuous
variables are indicated by the value ``'C'``. Binary variables should be specified
by ``'I'`` with lower bounds of 0 and upper bounds of 1.

A scalar is accepted for the ``sense``, ``b``, ``vartypes``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. 

The ``mixintprog`` function returns an instance of the type::
    
    type MixintprogSolution
        status
        objval
        sol
        attrs
    end

where ``status`` takes the same values as with ``linprog``.

If ``status`` does not indicate error or infeasiblity, the other members have the following values:

* ``objval`` -- optimal objective value
* ``sol`` -- primal solution vector
* ``attrs`` -- a dictionary that may contain other relevant attributes such as:

  - ``objbound`` -- Best known lower bound on the objective value

Analogous shortened and range-constraint versions are available as well.

We can solve a `binary knapsack problem <http://en.wikipedia.org/wiki/Knapsack_problem>`_

.. math::
    max\, &5x_1 + 3x_2 + 2x_3 + 7x_4 + 4x_5\\
    s.t.  &2x_1 + 8x_2 + 4x_3 + 2x_4 + 5x_5 \leq 10\\
          & (x_1, x_2, x_3, x_4, x_5) \in \{0,1\}^5

with the code::
    
    mixintprog(-[5.,3.,2.,7.,4.],[2. 8. 4. 2. 5.],'<',10,'I',0,1)

