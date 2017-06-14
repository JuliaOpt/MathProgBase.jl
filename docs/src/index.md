# MathProgBase

Solver-independent functions and low-level interfaces for Mathematical Programming

`MathProgBase.jl` <https://github.com/JuliaOpt/MathProgBase.jl> provides high-level one-shot functions for linear and mixed-integer programming, as well as a solver-independent low-level interface for implementing advanced techniques that require efficiently solving a sequence of linear programming problems.

To use MathProgBase, an external solver must be installed. See [Choosing Solvers](@ref).

```@contents
pages=[
    "linearprog.md",
    "mixedintprog.md",
    "quadprog.md",
    "solverinterface.md",
    "choosingsolver.md"
    ]
```
