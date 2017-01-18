----------------------
LinearQuadratic models
----------------------

This section describes the interface for ``LinearQuadratic`` models. This interface is for solvers which solve linear and quadratic programming problems and accept data as matrices which define the linear and quadratic components of the constraints and objective function. The concept of this interface similar to that of `OSI <https://projects.coin-or.org/Osi>`_, a C++ library which provides a generic virtual base class for interacting with linear programming solvers.

The interface allows for *hotstarts* by modifying some problem data and by adding constraints sequentially, if the solver supports them. We also provide an abstraction over *callbacks* for mixed-integer programming.


.. function:: loadproblem!(m::AbstractLinearQuadraticModel, filename::String)

    Loads problem data from the given file. Supported file types are solver-dependent.

.. function:: loadproblem!(m::AbstractLinearQuadraticModel, A, l, u, c, lb, ub, sense)

Loads the provided problem data to set up the linear programming problem:

.. math::
    \min_{x}\, &c^Tx\\
    s.t.     &lb \leq Ax \leq ub\\
             &l \leq x \leq u\\

``sense`` specifies the direction of the optimization problem, and must be either ``:Min`` or ``:Max``.

Both sparse and dense matrices are accepted for ``A``. ``Inf`` and ``-Inf`` indicate that
there is no corresponding upper or lower bound. Equal lower and upper bounds are used
to indicate equality constraints.

.. function:: writeproblem(m::AbstractLinearQuadraticModel, filename::String)

    Writes the current problem data to the given file. Supported file types are solver-dependent.



.. function:: getvarLB(m::AbstractLinearQuadraticModel)

    Returns a vector containing the lower bounds :math:`l` on the variables.

.. function:: setvarLB!(m::AbstractLinearQuadraticModel, l)

    Sets the lower bounds on the variables.

.. function:: getvarUB(m::AbstractLinearQuadraticModel)

    Returns a vector containing the upper bounds :math:`u` on the variables.

.. function:: setvarUB!(m::AbstractLinearQuadraticModel, u)

    Sets the upper bounds on the variables.

.. function:: getconstrLB(m::AbstractLinearQuadraticModel)

    Returns a vector containing the lower bounds :math:`lb` on the linear constraints.

.. function:: setconstrLB!(m::AbstractLinearQuadraticModel, lb)

    Sets the lower bounds on the linear constraints.

.. function:: getconstrUB(m::AbstractLinearQuadraticModel)

    Returns a vector containing the upper bounds :math:`ub` on the linear constraints.

.. function:: setconstrUB!(m::AbstractLinearQuadraticModel, ub)

    Sets the upper bounds on the linear constraints.

.. function:: getobj(m::AbstractLinearQuadraticModel)

    Returns a vector containing the linear objective coefficients :math:`c`.

.. function:: setobj!(m::AbstractLinearQuadraticModel, c)

    Sets the linear objective coefficients.

.. function:: getconstrmatrix(m::AbstractLinearQuadraticModel)

    Returns the full linear constraint matrix :math:`A`, typically as a
    ``SparseMatrixCSC``.

.. function:: addvar!(m::AbstractLinearQuadraticModel, constridx, constrcoef, l, u, objcoef)

    Adds a new variable to the model, with lower bound ``l`` (``-Inf`` if none),
    upper bound ``u`` (``Inf`` if none), and
    objective coefficient ``objcoef``. Constraint coefficients for this new variable
    are specified in a sparse format: the ``constrcoef`` vector contains the nonzero
    coefficients, and the ``constridx`` vector contains the indices of the corresponding
    linear constraints.

.. function:: addvar!(m::AbstractLinearQuadraticModel, l, u, objcoef)

    Adds a new variable to the model, with lower bound ``l`` (``-Inf`` if none),
    upper bound ``u`` (``Inf`` if none), and
    objective coefficient ``objcoef``. This is equivalent to calling the
    above method with empty arrays for the constraint coefficients.

.. function:: delvars!(m::AbstractLinearQuadraticModel, idxs)

    Removes the variables with indices in the vector idxs. The remaining variables in the model
    are renumbered. Let, for example, before deletion there be six variables with indices 
    1, 2, 3, 4, 5, 6, suppose we delete indices 1, 2, and 4, then there will be 3 variables remaining
    in the model, those which were previously numbered 3, 5, and 6. Now these 3 remaining variables
    respectly have indices 1, 2 and 3.

.. function:: addconstr!(m::AbstractLinearQuadraticModel, varidx, coef, lb, ub)

    Adds a new linear constraint to the model, with lower bound ``lb`` (``-Inf`` if none)
    and upper bound ``ub`` (``Inf`` if none). Coefficients for this new constraint
    are specified in a sparse format: the ``coef`` vector contains the nonzero
    coefficients, and the ``varidx`` vector contains the indices of the corresponding
    variables.

.. function:: delconstrs!(m::AbstractLinearQuadraticModel, idxs)

    Removes the constraints with indices in the vector idxs. The remaining constraints in the model
    are renumbered. Let, for example, before deletion there be six constraints with indices 
    1, 2, 3, 4, 5, 6, suppose we delete indices 1, 2, and 4, then there will be 3 constraints remaining
    in the model, those which were previously numbered 3, 5, and 6. Now these 3 remaining constraints 
    respectly have indices 1, 2 and 3.

.. function:: changecoeffs!(m::AbstractLinearQuadraticModel, cidxs, vidxs, val)

    Changes multiple coefficients in the A matrix. Coefficients to be changed are the ones with 
    constraint, variables and values indexed, respectively, in the vectors cidxs, vidxs and val.
    All the vectors must have the same size.

.. function:: numlinconstr(m::AbstractLinearQuadraticModel)

    Returns the number of linear constraints in the model.

.. function:: getconstrsolution(m::AbstractLinearQuadraticModel)

    Returns a vector containing the values of the linear constraints
    at the solution. This is the vector :math:`Ax`.

.. function:: getreducedcosts(m::AbstractLinearQuadraticModel)

    Returns the dual solution vector corresponding to the variable bounds,
    known as the reduced costs. Not available when integer variables are present.

.. function:: getconstrduals(m::AbstractLinearQuadraticModel)

    Returns the dual solution vector corresponding to the linear constraints.
    Not available when integer variables are present.

.. function:: getinfeasibilityray(m::AbstractLinearQuadraticModel)

    Returns a "Farkas" proof of infeasibility, i.e., an unbounded ray of the dual,
    for the linear constraints.
    Note that for some solvers, one must specify additional options for this
    ray to be computed.

.. function:: getbasis(m::AbstractLinearQuadraticModel)

    Returns the basis set for the optimal solution in the form ``(cbasis,rbasis)``,
    where both return values are vectors of symbols. The vector ``cbasis`` indexes
    the columns of the constraint matrix, while ``rbasis`` indexes the rows (values
    indicate whether the constraint is active at a lower/upper bound). The entries
    take value ``:Basic`` if the element is basic, ``:NonbasicAtLower`` if it is
    nonbasic at a lower bound, and ``:NonbasicAtUpper`` if it is nonbasic at upper
    bound. Other values may appear, taking solver-specific values. Note that this
    function may not work if the optimization algorithm is not able to provide
    basis information.

.. function:: getunboundedray(m::AbstractLinearQuadraticModel)

    Returns an unbounded ray of the problem, i.e., an objective-improving direction
    in which one may travel an infinite distance without violating any constraints.
    Note that for some solvers, one must specify additional options for this
    ray to be computed.

.. function:: getsimplexiter(m::AbstractLinearQuadraticModel)

    Returns the cumulative number of simplex iterations during the optimization process.
    In particular, for a MIP it returns the total simplex iterations for all nodes.

.. function:: getbarrieriter(m::AbstractLinearQuadraticModel)

    Returns the cumulative number of barrier iterations during the optimization process.

Integer Programming
^^^^^^^^^^^^^^^^^^^

.. function:: getnodecount(m::AbstractLinearQuadraticModel)

    Returns the total number of branch-and-bound nodes explored during the MIP optimization process.


.. function:: addsos1!(m::AbstractLinearQuadraticModel, idx, weight)

    Adds a special ordered set (SOS) constraint of type 1. Of the variables indexed by ``idx``, at most one can be nonzero. The ``weight`` argument induces the ordering of the variables; as such, they should be unique values. A typical SOS1 constraint might look like :math:`y=\sum_i w_i x_i`, where :math:`x_i \in \{0,1\}` are binary variables and the :math:`w_i` are weights. See `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for a description of SOS constraints and their potential uses.

.. function:: addsos2!(m::AbstractLinearQuadraticModel, idx, weight)

    Adds a special ordered set (SOS) constraint of type 2. Of the variables indexed by ``idx``, at most two can be nonzero, and if two are nonzero, they must be adjacent in the set. The ``weight`` argument induces the ordering of the variables; as such, they should be unique values. A common application for SOS2 constraints is modeling nonconvex piecewise linear functions; see `here <http://lpsolve.sourceforge.net/5.5/SOS.htm>`_ for details.

Quadratic Programming
^^^^^^^^^^^^^^^^^^^^^

.. function:: numquadconstr(m::AbstractLinearQuadraticModel)

    Returns the number of quadratic constraints in the model.

.. function:: setquadobj!(m::AbstractLinearQuadraticModel,Q)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. The argument ``Q`` must be either a symmetric positive semidefinite matrix or the upper triangular portion of a symmetric positive semidefinite matrix (when minimizing). Sparse (CSC) or dense representations are accepted.

.. function:: setquadobj!(m::AbstractLinearQuadraticModel,rowidx,colidx,quadval)

    Adds a quadratic term :math:`\frac{1}{2}x^TQx` to the objective, replacing any existing quadratic terms. Note the implicit :math:`\frac{1}{2}` scaling factor. Here the entries of :math:`Q` should be provided in sparse triplet form; e.g. entry indexed by ``k`` will fill ``quadval[k]`` in the ``(rowidx[k],colidx[k])`` entry of matrix ``Q``. Duplicate index sets ``(i,j)`` are accepted and will be summed together. Off-diagonal entries will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms. For example, ``setquadobj!(m, [1,1,2,2], [1,2,1,2], [3,1,1,1])`` and ``setquadobj!(m, [1,1,2], [1,2,2], [3,2,1])`` are both are valid descriptions for the matrix :math:`Q = \begin{pmatrix} 3 & 2 \\ 2 & 1 \end{pmatrix}`.

.. function:: setquadobjterms!(m::AbstractLinearQuadraticModel,rowidx,colidx,quadval)

    Provides an alternative "terms"-based interface to ``setquadobj!``. A list of quadratic terms is specified instead of the matrix ``Q``. For example, the objective :math:`x_1^2 + 2x_1x_2` is specified by ``setquadobjterms!(m,[1,1],[1,2],[1.0,2.0])``. Duplicate terms are summed together. Note: this method does not need to be implemented by solvers.

.. function:: addquadconstr!(m::AbstractLinearQuadraticModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)

    Adds the quadratic constraint :math:`s^Tx + \sum_{i,j} q_{i,j}x_ix_j \,\, sense \, rhs` to the model. The ``linearidx`` and ``linearval`` arrays specify the sparse vector ``s``. The quadratic terms are specified as in ``setquadobjterms!`` in the "terms" format. Sense must be ``'<'``, ``'>'``, or ``'='``. If supported by the solver, ``addquadconstr!`` may also be used to specify second-order cone (SOCP) and rotated second-order cone constraints. These should be of the form :math:`x^Tx -y^2 \le 0` or :math:`x^Tx -yz \le 0`, where :math:`y` and :math:`z` are restricted to be non-negative (in particular, :math:`Q` can have at most one off-diagonal term).

.. function:: getquadconstrsolution(m::AbstractLinearQuadraticModel)

    Returns a vector containing the values of the quadratic constraints
    at the solution.

.. function:: getquadconstrduals(m::AbstractLinearQuadraticModel)

    Returns the Lagrangian dual solution vector corresponding to the
    quadratic constraints. Some solvers do not compute these values by
    default. Not available when integer variables are present.

.. function:: getquadinfeasibilityray(m::AbstractLinearQuadraticModel)

    Returns a "Farkas" proof of infeasibility, i.e., an unbounded ray of the dual,
    for the quadratic constraints.
    Note that for some solvers, one must specify additional options for this
    ray to be computed.

.. function:: getquadconstrRHS(m::AbstractLinearQuadraticModel)

    Returns a vector containing the right-hand side values on the quadratic constraints.

.. function:: setquadconstrRHS!(m::AbstractLinearQuadraticModel, lb)

    Sets the right-hand side values on the quadratic constraints. If the constraint was provided in the special second-order conic format, the solver may reject changing the right-hand side from zero.
