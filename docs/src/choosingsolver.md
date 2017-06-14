# Choosing Solvers

Solvers and solver-specific parameters are specified by [`AbstractMathProgSolver`](@ref) objects, which are provided by particular solver packages. For example, the `Clp` package exports a `ClpSolver` object, which can be passed to [`linprog`](@ref) as follows::

```julia
    using Clp
    linprog([-1,0],[2 1],'<',1.5, ClpSolver())
```

Options are passed as keyword arguments, for example, `ClpSolver(LogLevel=1)`. See the [Clp](https://github.com/mlubin/Clp.jl), [Cbc](https://github.com/mlubin/Cbc.jl), [GLPKMathProgInterface](https://github.com/JuliaOpt/GLPKMathProgInterface.jl), and [Gurobi](https://github.com/JuliaOpt/Gurobi.jl) packages for more information.
