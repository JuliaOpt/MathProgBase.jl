.. _choosing-solvers:

----------------
Choosing solvers
----------------

Solvers and solver-specific parameters are specified by ``AbstractMathProgSolver`` objects, which are provided by particular solver packages. For example, the ``Clp`` package exports a ``ClpSolver`` object, which can be passed to ``linprog`` as follows::

    using Clp
    linprog([-1,0],[2 1],'<',1.5, ClpSolver())

Options are passed as keyword arguments, for example, ``ClpSolver(LogLevel=1)``. See the `Clp <https://github.com/mlubin/Clp.jl>`_, `Cbc <https://github.com/mlubin/Cbc.jl>`_, `GLPKMathProgInterface <https://github.com/JuliaOpt/GLPKMathProgInterface.jl>`_, and `Gurobi <https://github.com/JuliaOpt/Gurobi.jl>`_ packages for more information.
