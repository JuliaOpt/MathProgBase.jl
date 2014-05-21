------------------------
Semidefinite Programming
------------------------
MathProgBase supports a standardized interface for semidefinite programming. 

.. function:: addsdpvar!(m::AbstractMathProgModel, dim)
    
    Adds a semidefinite matrix variable to the model of dimension ``dim`` (:math:`S_+^{dim\times dim}`). Returns an integer index for referencing the variable in constraints or the objective.

.. function:: addsdpmatrix!(m::AbstractMathProgModel, mat)

    Adds a (symmetric) data matrix ``mat`` to the model. Can be used as coefficient in the objective or in a constraint. Returns an integer index for referencing the matrix in constraints or the objective.

.. function:: addsdpconstr!(m::AbstractMathProgModel, matvaridx, matcoefidx, scalidx, scalcoef, lb, ub)

    Adds a semidefinite constraint in "primal" form with optional scalar variable terms, i.e. :math:`ub \leq \sum_{i} A_i \cdot X_i + \sum_{j} c_j y_j \leq ub`. ``matvaridx`` (resp. ``scalidx``) are the indices of the matrix (resp. scalar) variables ``X_i`` (resp. ``y_j``) in the model, while ``matcoefidx`` are the indices of the data matrices ``A_i`` in the model. ``scalcoef`` are the coefficients for the scalar variable terms ``c_j``. Returns an integer index for querying dual values after optimization via ``getsdpdual``.

.. function:: setsdpobj!(m::AbstractMathProgModel, matvaridx, matcoefidx)

    Sets semidefinite terms in the objective. Data and variables are passed via ``matvaridx`` and ``matcoefidx`` as described in ``addsdpconstr!``. Note that scalar terms can be added in conjunction via ``setobj!``.

.. function:: getsdpsolution(m::AbstractMathProgModel, idx)

    Grabs the solution value for semidefinite matrix ``idx`` after optimization.

.. function:: getsdpdual(m::AbstractMathProgModel)
    
    Grabs dual solutions for semidefinite constraints. Returns a vector with entries corresponding to the indices returned by ``addsdpconstr!``.
