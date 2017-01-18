----------------------
SolverInterface module
----------------------

The high-level functions ``linprog``, ``mixintprog``, and ``quadprog`` are
are written on top of a solver-independent low-level interface defined in the
``SolverInterface`` module. This is also the abstraction layer that modeling languages
`JuMP <https://github.com/JuliaOpt/JuMP.jl>`_ and
`Convex.jl <https://github.com/JuliaOpt/Convex.jl>`_ use to communicate with solvers.

The intent of the ``SolverInterface`` module is to provide an abstraction over
the low-level interfaces which are common among solvers. By grouping solvers
into a small number of categories, we can more easily implement automated
transformations between different representations of problems,
providing a mapping between the user's representation of the problem
and the data structures which the solver takes natively as input.

Solvers are divided into three categories (some solvers belong to multiple categories):

- ``LinearQuadratic``: These are solvers which solve linear and quadratic programming problems and accept data as matrices which define the linear and quadratic components of the constraints and objective function. Some of these solvers support warm-starting when certain data are changed or when constraints are added (e.g., by reusing the previous optimal basis from the simplex method). These include classical mixed-integer linear programming solvers with support for callbacks. Some solvers support adding second-order conic constraints through specially formatted quadratic constraints. Examples of solvers in this category are `Cbc <https://github.com/JuliaOpt/Cbc.jl>`_, `Clp <https://github.com/JuliaOpt/Clp.jl>`_, `CPLEX <https://github.com/JuliaOpt/CPLEX.jl>`_, `GLPK <https://github.com/JuliaOpt/GLPK.jl>`_, `Gurobi <https://github.com/JuliaOpt/Gurobi.jl>`_, and `Mosek <https://github.com/JuliaOpt/Mosek.jl>`_.
- ``Conic``: These are solvers which solve conic programming problems; these are problems with linear objectives where affine functions of decision variables are restricted to fall within certain convex cones (e.g., the second-order cone, positive semidefinite cone, and exponential cone). The input format for these solves are the matrices and vectors defining the affine functions and the list of cones. Examples of solvers in this category are `ECOS <https://github.com/JuliaOpt/ECOS.jl>`_, `Mosek <https://github.com/JuliaOpt/Mosek.jl>`_, and `SCS <https://github.com/JuliaOpt/SCS.jl>`_.
- ``Nonlinear``: These are traditional derivative-based nonlinear solvers which require derivative evaluation oracles, as well as solvers which require access to the algebraic representation of a problem (when available). Examples of solvers in this category are `AmplNLWriter <https://github.com/JuliaOpt/AmplNLWriter.jl>`_, `CoinOptServices <https://github.com/JuliaOpt/CoinOptServices.jl>`_, `Ipopt <https://github.com/JuliaOpt/Ipopt.jl>`_, `KNITRO <https://github.com/JuliaOpt/KNITRO.jl>`_, `Mosek <https://github.com/JuliaOpt/Mosek.jl>`_, and `NLopt <https://github.com/JuliaOpt/NLopt.jl>`_.


The ``SolverInterface`` module makes a distinction between a "solver" and a "model". A solver is a lightweight object used for selecting solvers and parameters. It does not store any instance data. Solver objects are used to create new instances of model objects. Model objects should be seen as a solver's in-memory representation of a particular problem.

Solver objects should inherit from the ``AbstractMathProgSolver`` abstract type. Solvers should implement at least one of the following methods:

.. function:: LinearQuadraticModel(s::AbstractMathProgSolver)

    Returns an instance of an ``AbstractLinearQuadraticModel`` using the given solver.

.. function:: ConicModel(s::AbstractMathProgSolver)

    Returns an instance of an ``AbstractConicModel`` using the given solver.

.. function:: NonlinearModel(s::AbstractMathProgSolver)

    Returns an instance of an ``AbstractNonlinearModel`` using the given solver.

The interfaces for these three classes of models are defined in the subsequent sections.

All abstract model types inherit from the abstract type ``AbstractMathProgModel``. The following methods are common to all model types:

.. function:: getsolution(m::AbstractMathProgModel)

    Returns the solution vector found by the solver.

.. function:: getobjval(m::AbstractMathProgModel)

    Returns the objective value of the solution found by the solver. In particular, this may be
    the objective value for the best feasible solution if optimality is not attained or proven.

.. function:: optimize!(m::AbstractMathProgModel)

    Solves the optimization problem.

.. function:: status(m::AbstractMathProgModel)

    Returns the termination status after solving (i.e., the reason why the solver terminated). Possible values include ``:Optimal``,
    ``:Infeasible``, ``:Unbounded``, ``:DualityFailure``, ``:UserLimit`` (iteration limit or timeout), and ``:Error``.
    Solvers may return other statuses, for example, when presolve indicates that the model is
    either infeasible or unbounded, but did not determine which.

.. function:: getobjbound(m::AbstractMathProgModel)

    Returns the best known bound on the optimal objective value.
    This is used, for example, when a branch-and-bound method
    is stopped before finishing.

.. function:: getobjgap(m::AbstractMathProgModel)

    Returns the final relative optimality gap as optimization terminated. That is, it returns
    :math:`\frac{|b-f|}{|f|}`, where :math:`b` is the best bound and :math:`f` is the best
    feasible objective value.

.. function:: getrawsolver(m::AbstractMathProgModel)

    Returns an object that may be used to access a solver-specific API for this model.

.. function:: getsolvetime(m::AbstractMathProgModel)

    Returns the total elapsed solution time as reported by the solver.

.. function:: setsense!(m::AbstractMathProgModel, sense)

    Sets the optimization sense of the model. Accepted values are ``:Min`` and ``:Max``.

.. function:: getsense(m::AbstractMathProgModel)

    Returns the optimization sense of the model.

.. function:: numvar(m::AbstractMathProgModel)

    Returns the number of variables in the model.

.. function:: numconstr(m::AbstractMathProgModel)

    Returns the total number of constraints in the model.

.. function:: freemodel!(m::AbstractMathProgModel)

    Release any resources and memory used by the model. Note that the
    Julia garbage collector takes care of this automatically, but
    automatic collection cannot always be forced. This method is useful for more
    precise control of resources, especially in the case of commercial solvers
    with licensing restrictions on the number of concurrent runs.
    Users must discard the model object after this method is invoked.

.. function:: copy(m::AbstractMathProgModel)

    Copies the model m by replicating the solver internal model, callbacks are not copied.

.. function:: setvartype!(m::AbstractMathProgModel, v::Vector{Symbol})

    Sets the types of the variables to those indicated by the vector ``v``. Valid
    types are ``:Int`` for integer, ``:Cont`` for continuous, ``:Bin`` for binary,
    ``:SemiCont`` for `semicontinuous <http://orinanobworld.blogspot.com/2011/03/semicontinuous-variables.html>`_, and ``:SemiInt`` for `semi-integer <http://www.gams.com/mccarl/mccarlhtml/semi-integer_variables.htm>`_.

.. function:: getvartype(m::AbstractMathProgModel)

    Returns a vector indicating the types of each variable, with values described above.

It is the philosophy of MathProgBase to not abstract over most solver parameters, because very few are universal, and we do not want to make it difficult for users who already familiar with the parameters of their solver of choice. However, in certain situations an abstraction over parameters is needed, for example, in meta-solvers which must adjust parameters of a subsolver inside a loop. Parameters set using these methods should override any parameters provided in a solver-dependent way. Solvers/models may chose to implement the following method:

.. function:: setparameters!(m::Union{AbstractMathProgSolver,AbstractMathProgModel}; kwargs...)

    Sets solver-independent parameters via keyword arguments. Curent valid parameters are ``TimeLimit`` and ``Silent``.
    ``TimeLimit`` (``Float64``): If the solve is not completed to optimality tolerances within the given number of seconds, the solver
    should return immediately with status ``:UserLimit``.
    ``Silent`` (``Bool``): If set to true, the solver should disable any output to the console. If set to false, the parameter
    has no effect.

If these parameter-setting methods are called on an ``AbstractMathProgSolver``, then they should apply to all new models created from the solver (but not existing models). If they are called on an ``AbstractMathProgModel``, they should apply to that model only. Unrecognized parameters (those not listed above) should be ignored or trigger a warning message.

.. function:: setwarmstart!(m::AbstractMathProgModel, v)

    Provide an initial solution ``v`` to the solver, as supported. To leave values undefined, set them
    to ``NaN``. MIP solvers should ignore provided solutions that are infeasible or
    cannot be completed to a feasible solution. Nonlinear solvers may use provided
    solutions as starting points even if infeasible.

