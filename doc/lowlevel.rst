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

.. function:: getconstrmatrix(m::AbstractMathProgModel)

    Returns the full constraint matrix :math:`A`, typically as a 
    ``SparseMatrixCSC``.

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

.. function:: getbasis(m::AbstractMathProgModel)

    Returns the basis set for the optimal solution in the form ``(cbasis,rbasis)``, 
    where both return values are vectors of symbols. The vector ``cbasis`` indexes 
    the columns of the constraint matrix, while ``rbasis`` indexes the rows (values 
    indicate whether the constraint is active at a lower/upper bound). The entries
    take value ``:Basic`` if the element is basic, ``:NonbasicAtLower`` if it is 
    nonbasic at a lower bound, and ``:NonbasicAtUpper`` if it is nonbasic at upper 
    bound. Other values may appear, taking solver-specific values. Note that this 
    function may not work if the optimization algorithm is not able to provide 
    basis information.

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

.. function:: addsos1!(m::AbstractMathProgModel, idx, weight)
    
    Adds a special ordered set (SOS) constraint of type 1. Of the variables indexed by ``idx``, at most one can be nonzero. The ``weight`` argument induces the ordering of the variables; as such, they should be unique values. A typical SOS1 constraint might look like :math:`y=\sum_i w_i x_i`, where :math:`x_i \in \{0,1\}` are binary variables and the :math:`w_i` are weights. See `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for a description of SOS constraints and their potential uses.

.. function:: addsos2!(m::AbstractMathProgModel, idx, weight)
    
    Adds a special ordered set (SOS) constraint of type 2. Of the variables indexed by ``idx``, at most two can be nonzero, and if two are nonzero, they must be adjacent in the set. The ``weight`` argument induces the ordering of the variables; as such, they should be unique values. A common application for SOS2 constraints is modeling nonconvex piecewise linear functions; see `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for details.

.. function:: setquadobj!(m::AbstractMathProgModel,Q)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. The argument ``Q`` must be either a symmetric positive semidefinite matrix or the upper triangular portion of a symmetric positive semidefinite matrix (when minimizing). Sparse (CSC) or dense representations are accepted.

.. function:: setquadobj!(m::AbstractMathProgModel,rowidx,colidx,quadval)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. The matrix :math:`Q` must be symmetric positive semidefinite (when minimizing). Here the entries of :math:`Q` should be provided in sparse triplet form; e.g. entry indexed by ``k`` will fill ``quadval[k]`` in the ``(rowidx[k],colidx[k])`` entry of matrix ``Q``. Duplicate index sets ``(i,j)`` are accepted and will be summed together. Off-diagonal entries will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms. For example, ``setquadobj!(m, [1,1,2,2], [1,2,1,2], [3,1,1,1])`` and ``setquadobj!(m, [1,1,2], [1,2,2], [3,2,1])`` are both are valid descriptions for the matrix :math:`Q = \begin{pmatrix} 3 & 2 \\ 2 & 1 \end{pmatrix}`.

.. function:: setquadobjterms!(m::AbstractMathProgModel,rowidx,colidx,quadval)

    Provides an alternative "terms"-based interface to ``setquadobj!``. A list of quadratic terms is specified instead of the matrix ``Q``. For example, the objective :math:`x_1^2 + 2x_1x_2` is specified by ``setquadobjterms!(m,[1,1],[1,2],[1.0,2.0])``. Duplicate terms are summed together. Note: this method does not need to be implemented by solvers.

.. function:: addquadconstr!(m::AbstractMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)

    Adds the quadratic constraint :math:`s^Tx + \sum_{i,j} q_{i,j}x_ix_j \,\, sense \, rhs` to the model. The ``linearidx`` and ``linearval`` arrays specify the sparse vector ``s``. The quadratic terms are specified as in ``setquadobjterms!`` in the "terms" format. Sense must be ``'<'`` or ``'>'``, and :math:`Q` must be positive semidefinite or negative semidefinite, respectively. If supported by the solver, ``addquadconstr!`` may also be used to specify second-order cone (SOCP) and rotated second-order cone constraints. These should be of the form :math:`x^Tx -y^2 \le 0` or :math:`x^Tx -yz \le 0`, where :math:`y` and :math:`z` are restricted to be non-negative (in particular, :math:`Q` can have at most one off-diagonal term).
