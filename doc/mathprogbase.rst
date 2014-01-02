===============
MathProgBase.jl
===============

.. module:: MathProgBase
   :synopsis: Solver-independent functions and low-level interfaces for Mathematical Programming

`MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl>`_ provides high-level
one-shot functions for linear and mixed-integer programming,
as well as a solver-independent low-level interface for implementing advanced techniques
that require efficiently solving a sequence of linear programming problems.

To use MathProgBase, an external solver must be installed. See :ref:`choosing solvers <choosing-solvers>`.

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
*    ``A`` is the constraint matrix, with rows :math:`a_i`
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` is an *optional* parameter specifying the desired solver, see :ref:`choosing solvers <choosing-solvers>`. If this parameter is not provided, the default solver is used.
 
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

    using MathProgBase
    
    sol = linprog([-1,0],[2 1],'<',1.5)
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
*    ``solver`` is an *optional* parameter specifying the desired solver, see :ref:`choosing solvers <choosing-solvers>`. If this parameter is not provided, the default solver is used.
 
A scalar is accepted for the ``l``, ``u``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

A shortened version is defined as::

    linprog(c, A, lb, ub, solver) = linprog(c, A, lb, ub, 0, Inf, solver)


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

---------------------
Quadratic Programming
---------------------

.. function:: quadprog(c, Q, A, sense, b, l, u, solver)

Solves the quadratic programming problem:

.. math::
    \min_{x}\, &\frac{1}{2}x^TQx + c^Tx\\
    s.t.     &a_i^Tx \text{ sense}_i \, b_i \forall\,\, i\\
             &l \leq x \leq u\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``Q`` is the Hessian matrix of the objective
*    ``A`` is the constraint matrix, with rows :math:`a_i`
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` is an *optional* parameter specifying the desired solver, see :ref:`choosing solvers <choosing-solvers>`. If this parameter is not provided, the default solver is used.
 
A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. note::
    Quadratic programming solvers extensively exploit the sparsity of the Hessian matrix ``Q`` and the constraint matrix ``A``. While both dense and sparse matrices are accepted, for large-scale problems sparse matrices should be provided if permitted by the problem structure.  

The ``quadprog`` function returns an instance of the type::
    
    type QuadprogSolution
        status
        objval
        sol
        attrs
    end

where ``status`` is a termination status symbol, one of ``:Optimal``, ``:Infeasible``, ``:Unbounded``, ``:UserLimit`` (iteration limit or timeout), ``:Error`` (and maybe others).

If ``status`` is ``:Optimal``, the other members have the following values:

* ``objval`` -- optimal objective value
* ``sol`` -- primal solution vector
* ``attrs`` -- a dictionary that may contain other relevant attributes (not currently used).

Analogous shortened and range-constraint versions are available as well.

We can solve the three-dimensional QP (see ``test/quadprog.jl``):

.. math::
    \min_{x,y,z}\, &x^2+y^2+z^2+xy+yz\\
    s.t.         &x + 2y + 3z \leq 4\\
                 &x + y \leq 1

by::

    using MathProgBase
    
    sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],[4., 1.],Inf,-Inf,Inf)
    if sol.status == :Optimal
        println("Optimal objective value is $(sol.objval)")
        println("Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2]), $(sol.sol[3])]")
    else
        println("Error: solution status $(sol.status)")
    end

-------------------
Low-level interface
-------------------

The ``linprog`` and ``mixintprog`` functions are written on top of a solver-independent low-level interface called ``MathProgSolverInterface``, which individual solvers implement. The concept is similar to that of `OSI <https://projects.coin-or.org/Osi>`_, a C++ library which provides a generic virtual base class for interacting with linear programming solvers. Julia, however, does not quite have virtual classes or interfaces. Instead, multiple dispatch is used with abstract types. The API is designed to support problem modification as needed to solve a sequence of linear programming problems efficiently; linear programming solvers are expected to hot-start the solution process after modifications such as additional constraints or variables. For mixed-integer programming, hot-starting is usually impractical.

The ``MathProgSolverInterface`` exports two abstract types: ``AbstractMathProgModel``, which represents an instance of an optimization problem, and ``AbstractMathProgSolver``, which represents a solver (with particular solution options), from which an ``AbstractMathProgModel`` is generated.

.. function:: model(s::AbstractMathProgSolver)

    Returns an instance of an ``AbstractMathProgModel`` using the given solver.

.. function:: loadproblem!(m::AbstractMathProgModel, filename::String)

    Loads problem data from the given file. Supported file types are solver-dependent.

.. function:: loadproblem!(m::AbstractMathProgModel, A, l, u, c, lb, ub, sense)

Loads the provided problem data to set up the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &lb \leq Ax \leq ub\\
             &l \leq x \leq u\\

``sense`` specifies the direction of the optimization problem, and must be either ``:Min`` or ``:Max``.

Both sparse and dense matrices are accepted for ``A``. ``Inf`` and ``-Inf`` indicate that 
there is no corresponding upper or lower bound. Equal lower and upper bounds are used
to indicate equality constraints.

.. function:: writeproblem(m::AbstractMathProgModel, filename::String)
    
    Writes the current problem data to the given file. Supported file types are solver-dependent.



.. function:: getvarLB(m::AbstractMathProgModel)
   
    Returns a vector containing the lower bounds :math:`l` on the variables.

.. function:: setvarLB!(m::AbstractMathProgModel, l)
   
    Sets the lower bounds on the variables.

.. function:: getvarUB(m::AbstractMathProgModel)
   
    Returns a vector containing the upper bounds :math:`u` on the variables.

.. function:: setvarUB!(m::AbstractMathProgModel, u)
   
    Sets the upper bounds on the variables.



.. function:: getconstrLB(m::AbstractMathProgModel)
   
    Returns a vector containing the lower bounds :math:`lb` on the constraints.

.. function:: setconstrLB!(m::AbstractMathProgModel, lb)
   
    Sets the lower bounds on the constraints.

.. function:: getconstrUB(m::AbstractMathProgModel)
   
    Returns a vector containing the upper bounds :math:`ub` on the constraints.

.. function:: setconstrUB!(m::AbstractMathProgModel, ub)
   
    Sets the upper bounds on the constraints.

.. function:: getobj(m::AbstractMathProgModel)
   
    Returns a vector containing the objective coefficients :math:`c`.

.. function:: setobj!(m::AbstractMathProgModel, c)
   
    Sets the objective coefficients.



.. function:: addvar!(m::AbstractMathProgModel, constridx, constrcoef, l, u, objcoef)

    Adds a new variable to the model, with lower bound ``l`` (``-Inf`` if none), 
    upper bound ``u`` (``Inf`` if none), and
    objective coefficient ``objcoef``. Constraint coefficients for this new variable
    are specified in a sparse format: the ``constrcoef`` vector contains the nonzero
    coefficients, and the ``constridx`` vector contains the indices of the corresponding
    constraints.

.. function:: addvar!(m::AbstractMathProgModel, l, u, objcoef)

    Adds a new variable to the model, with lower bound ``l`` (``-Inf`` if none), 
    upper bound ``u`` (``Inf`` if none), and
    objective coefficient ``objcoef``. This is equivalent to calling the 
    above method with empty arrays for the constraint coefficients.
    

.. function:: addconstr!(m::AbstractMathProgModel, varidx, coef, lb, ub)

    Adds a new constraint to the model, with lower bound ``lb`` (``-Inf`` if none)
    and upper bound ``ub`` (``Inf`` if none). Coefficients for this new constraint
    are specified in a sparse format: the ``coef`` vector contains the nonzero
    coefficients, and the ``varidx`` vector contains the indices of the corresponding
    variables.



.. function:: updatemodel!(m::AbstractMathProgModel)

    Commits recent changes to the model. Only required by some solvers (e.g. Gurobi).

.. function:: setsense!(m::AbstractMathProgModel, sense)

    Sets the optimization sense of the model. Accepted values are ``:Min`` and ``:Max``.

.. function:: getsense(m::AbstractMathProgModel)

    Returns the optimization sense of the model.

.. function:: numvar(m::AbstractMathProgModel)

    Returns the number of variables in the model.

.. function:: numconstr(m::AbstractMathProgModel)

    Returns the number of constraints in the model.

.. function:: optimize!(m::AbstractMathProgModel)

    Solves the optimization problem.

.. function:: status(m::AbstractMathProgModel)

    Returns the termination status after solving. Possible values include ``:Optimal``,
    ``:Infeasible``, ``:Unbounded``, ``:UserLimit`` (iteration limit or timeout), and ``:Error``.
    Solvers may return other statuses, for example, when presolve indicates that the model is
    either infeasible or unbounded, but did not determine which.

.. function:: getobjval(m::AbstractMathProgModel)

    Returns the objective value of the solution found by the solver.

.. function:: getobjbound(m::AbstractMathProgModel)

    Returns the best known bound on the optimal objective value.
    This is used, for example, when a branch-and-bound method
    is stopped before finishing.

.. function:: getsolution(m::AbstractMathProgModel)

    Returns the solution vector found by the solver.

.. function:: getconstrsolution(m::AbstractMathProgModel)

    Returns a vector containing the values of the constraints
    at the solution. This is the vector :math:`Ax`.

.. function:: getreducedcosts(m::AbstractMathProgModel)

    Returns the dual solution vector corresponding to the variable bounds,
    known as the reduced costs. Not available when integer variables are present.

.. function:: getconstrduals(m::AbstractMathProgModel)

    Returns the dual solution vector corresponding to the constraints.
    Not available when integer variables are present.

.. function:: getinfeasibilityray(m::AbstractMathProgModel)

    Returns a "Farkas" proof of infeasibility, i.e., an unbounded ray of the dual. 
    Note that for some solvers, one must specify additional options for this
    ray to be computed.

.. function:: getunboundedray(m::AbstractMathProgModel)

    Returns an unbounded ray of the problem, i.e., an objective-improving direction 
    in which one may travel an infinite distance without violating any constraints.
    Note that for some solvers, one must specify additional options for this
    ray to be computed.

.. function:: getrawsolver(m::AbstractMathProgModel)

    Returns an object that may be used to access a solver-specific API for this model.

.. function:: setvartype!(m::AbstractMathProgModel, v::Vector{Char})

    Sets the types of the variables to those indicated by the vector ``v``. Valid
    types are ``'I'`` for integer and ``'C'`` for continuous. Binary variables
    should be indicated by ``'I'`` with lower bound 0 and upper bound 1.

.. function:: getvartype(m::AbstractMathProgModel)

    Returns a vector indicating the types of each variable, with values described above.

.. function:: setwarmstart!(m::AbstractMathProgModel, v)

    Provide an initial solution ``v`` to the MIP solver. To leave values undefined, set them
    to ``NaN``.

.. function:: setquadobj!(m::AbstractMathProgModel,Q)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. The argument ``Q`` must be either a symmetric positive semidefinite matrix or the upper triangular portion of a symmetric positive semidefinite matrix (when minimizing). Sparse (CSC) or dense representations are accepted.

.. function:: setquadobj!(m::AbstractMathProgModel,rowidx,colidx,quadval)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. The matrix :math:`Q` must be symmetric positive semidefinite (when minimizing). Here the entries of :math:`Q` should be provided in sparse triplet form; e.g. entry indexed by ``k`` will fill ``quadval[k]`` in the ``(rowidx[k],colidx[k])`` entry of matrix ``Q``. Duplicate index sets ``(i,j)`` are accepted and will be summed together. Off-diagonal entries will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms. For example, ``setquadobj!(m, [1,1,2,2], [1,2,1,2], [3,1,1,1])`` and ``setquadobj!(m, [1,1,2], [1,2,2], [3,2,1])`` are both are valid descriptions for the matrix :math:`Q = \begin{pmatrix} 3 & 2 \\ 2 & 1 \end{pmatrix}`.

.. function:: setquadobjterms!(m::AbstractMathProgModel,rowidx,colidx,quadval)

    Provides an alternative "terms"-based interface to ``setquadobj!``. A list of quadratic terms is specified instead of the matrix ``Q``. For example, the objective :math:`x_1^2 + 2x_1x_2` is specified by ``setquadobjterms!(m,[1,1],[1,2],[1.0,2.0])``. Duplicate terms are summed together. Note: this method does not need to be implemented by solvers.

.. function:: addquadconstr!(m::AbstractMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)

    Adds the quadratic constraint :math:`s^Tx + \frac{1}{2}x^TQx \,\, sense \, rhs` to the model. The ``linearidx`` and ``linearval`` arrays specify the sparse vector ``s``. The quadratic term is specified as in ``setquadobj!``. Sense must be ``'<'`` or ``'>'``, and :math:`Q` must be positive semidefinite or negative semidefinite, respectively. If supported by the solver, ``addquadconstr!`` may also be used to specify second-order cone (SOCP) and rotated second-order cone constraints. These should be of the form :math:`x^Tx -y^2 \le 0` or :math:`x^Tx -yz \le 0`, where :math:`y` and :math:`z` are restricted to be non-negative (in particular, :math:`Q` can have at most one off-diagonal term).


.. _choosing-solvers:

----------------
Choosing solvers
----------------

Solvers and solver-specific parameters are specified by ``AbstractMathProgSolver`` objects, which are provided by particular solver packages. For example, the ``Clp`` package exports a ``ClpSolver`` object, which can be passed to ``linprog`` as follows::

    using Clp
    linprog([-1,0],[2 1],'<',1.5, ClpSolver())

Options are passed as keyword arguments, for example, ``ClpSolver(LogLevel=1)``. See the `Clp <https://github.com/mlubin/Clp.jl>`_, `Cbc <https://github.com/mlubin/Cbc.jl>`_, `GLPKMathProgInterface <https://github.com/JuliaOpt/GLPKMathProgInterface.jl>`_, and `Gurobi <https://github.com/JuliaOpt/Gurobi.jl>`_ packages for more information.

If no solver is specified, a default is chosen. See ``src/defaultsolvers.jl`` for the list of default solvers.

-------------
MIP Callbacks
-------------
MathProgBase supports a standardized and abstracted way to implement common MIP callbacks on the model. Currently there is support for adding:

*    Lazy constraints (only added to model if violated by integer-feasible solution)
*    Cut callbacks (only cuts off non-integer feasible solutions)
*    Heuristic callbacks (proposes heuristically constructed integer-feasible solutions at MIP nodes)

A more detailed description of the three types of supported callbacks can be found in the JuMP documentation `here <https://jump.readthedocs.org/en/latest/jump.html#solver-callbacks>`_.

The ``MathProgSolverInterface`` exports an abstract type ``MathProgCallbackData`` which represents the solver-specific data needed to implement the callback.

.. function:: setlazycallback!(m::AbstractMathProgModel,f)

   Adds lazy constraint callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: setcutcallback!(m::AbstractMathProgModel,f)

   Adds cut callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: setheuristiccallback!(m::AbstractMathProgModel,f)

   Adds heuristic callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: cbgetmipsolution(d::MathProgCallbackData[, output])

   Grabs current best integer-feasible solution to the model. The optional second argument specifies an output vector.
   
.. function:: cbgetlpsolution(d::MathProgCallbackData[, output])

   Grabs current best linear relaxation solution to the model. The optional second argument specifies an output vector.
   
.. function:: cbgetobj(d::MathProgCallbackData)

   Grabs objective value for current best integer-feasible solution.
   
.. function:: cbgetbestbound(d::MathProgCallbackData) 

   Grabs best bound for objective function found so far (lower bound when minimizing, upper bound when maximizing).
   
.. function:: cbgetexplorednodes(d::MathProgCallbackData)

   Returns number of nodes that have been explored so far in the solve process.

.. function:: cbgetstate(d::MathProgCallbackData)

   Returns current location in solve process: ``:MIPNode`` if at node in branch-and-cut tree, ``:MIPSol`` at an integer-feasible solution, and ``:Other`` otherwise.

.. function:: cbaddsolution!(d::MathProgCallbackData,x)

   Adds feasible solution ``x`` to model.
   
.. function:: cbaddcut!(d::MathProgCallbackData,varidx,varcoef,sense,rhs) 

   Adds cut to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.
   
.. function:: cbaddlazy!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds lazy constraint to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.
   
