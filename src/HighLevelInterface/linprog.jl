type LinprogSolution
    status
    objval
    sol
    attrs
end

InputVector{T<:Union{Real,Char}} = Union{Vector{T},Real,Char}

function expandvec(x,len::Integer)
    if isa(x,Vector)
        if length(x) != len
            error("Input size mismatch. Expected vector of length $len but got $(length(x))")
        end
        return x
    else
        return fill(x,len)
    end
end

function no_lp_solver()
    error("""
          A linear programming solver must be specified as the final argument.

          You can use any solver listed on the solvers table of http://www.juliaopt.org/ that has a checkmark in the Linear/Quadratic column.

          A (free) example is Clp.jl. Once Clp is installed and imported via "using Clp", you can specify ClpSolver() as the solver.
          Solver options are specified by using keyword arguments to ClpSolver().
          """)
end

"""
    buildlp(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

    Builds the linear programming problem as defined in [`linprog`](@ref) and accepts the following arguments:

*    ``c`` is the objective vector, always in the sense of minimization

*    ``A`` is the constraint matrix

*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``

*    ``b`` is the right-hand side vector

*    ``l`` is the vector of lower bounds on the variables

*    ``u`` is the vector of upper bounds on the variables, and

*    ``solver`` specifies the desired solver, see [Choosing Solvers](@ref).

A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

This variant of [`buildlp`](@ref) allows to specify two-sided linear constraints (also known as range constraints) similar to [`linprog`](@ref), and accepts the following arguments:

    buildlp(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

*    ``c`` is the objective vector, always in the sense of minimization

*    ``A`` is the constraint matrix

*    ``rowlb`` is the vector of row lower bounds

*    ``rowub`` is the vector of row upper bounds

*    ``lb`` is the vector of lower bounds on the variables

*    ``ub`` is the vector of upper bounds on the variables, and

*    ``solver`` specifies the desired solver, see [Choosing Solvers](@ref).

A scalar is accepted for the ``l``, ``u``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

The [`buildlp`](@ref) function returns an [`AbstractLinearQuadraticModel`](@ref) that can be input to [`solvelp`](@ref) in order to obtain a solution.
"""

function buildlp(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)
    m = LinearQuadraticModel(solver)
    nrow,ncol = size(A)

    c = expandvec(c, ncol)
    rowlbtmp = expandvec(rowlb, nrow)
    rowubtmp = expandvec(rowub, nrow)
    lb = expandvec(lb, ncol)
    ub = expandvec(ub, ncol)

    # rowlb is allowed to be vector of senses
    if eltype(rowlbtmp) == Char
        realtype = eltype(rowubtmp)
        warn_no_inf(realtype)
        sense = rowlbtmp
        rhs = rowubtmp
        @assert realtype <: Real

        rowlb = Array{realtype}(nrow)
        rowub = Array{realtype}(nrow)
        for i in 1:nrow
            if sense[i] == '<'
                rowlb[i] = typemin(realtype)
                rowub[i] = rhs[i]
            elseif sense[i] == '>'
                rowlb[i] = rhs[i]
                rowub[i] = typemax(realtype)
            elseif sense[i] == '='
                rowlb[i] = rhs[i]
                rowub[i] = rhs[i]
            else
                error("Unrecognized sense '$(sense[i])'")
            end
        end
    else
        rowlb = rowlbtmp
        rowub = rowubtmp
    end

    loadproblem!(m, A, lb, ub, c, rowlb, rowub, :Min)

    return m
end

"""
Solves the linear programming problem as defined in [`linprog`](@ref) and accepts the following argument:

*    `m` is an [`AbstractLinearQuadraticModel`](@ref) (e.g., as returned by [`buildlp`](@ref)).

The [`solvelp`](@ref) function returns an instance of the type::

```julia
type LinprogSolution
    status
    objval
    sol
    attrs
end
```

"""
function solvelp(m)
    optimize!(m)
    stat = status(m)
    if stat == :Optimal
        attrs = Dict()
        attrs[:redcost] = getreducedcosts(m)
        attrs[:lambda] = getconstrduals(m)
        return LinprogSolution(stat, getobjval(m), getsolution(m), attrs)
    elseif stat == :Infeasible
        attrs = Dict()
        try
            attrs[:infeasibilityray] = getinfeasibilityray(m)
        catch
            warn("Problem is infeasible, but infeasibility ray (\"Farkas proof\") is unavailable; check that the proper solver options are set.")
        end
        return LinprogSolution(stat, nothing, [], attrs)
    elseif stat == :Unbounded
        attrs = Dict()
        try
            attrs[:unboundedray] = getunboundedray(m)
        catch
            warn("Problem is unbounded, but unbounded ray is unavailable; check that the proper solver options are set.")
        end
        return LinprogSolution(stat, nothing, [], attrs)
    else
        return LinprogSolution(stat, nothing, [], Dict())
    end
end
"""

    linprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

This function allows one to specify two-sided linear constraints (also known as range constraints) to solve the linear programming problem:

```math
\\begin{align*}
\\text{min}_{x} \\quad & c^{T} x \\\\
\\text{s.t.}    \\quad & rowlb \\leq A^{T} x \\leq rowub \\\\
                \\quad & l \\leq x \\leq u \\\\

\\end{align*}
```

where:

*    ``c`` is the objective vector, always in the sense of minimization

*    ``A`` is the constraint matrix

*    ``rowlb`` is the vector of row lower bounds

*    ``rowub`` is the vector of row upper bounds

*    ``lb`` is the vector of lower bounds on the variables

*    ``ub`` is the vector of upper bounds on the variables, and

*    ``solver`` specifies the desired solver, see [Choosing Solvers](@ref).

A scalar is accepted for the ``l``, ``u``, ``rowlb``, and ``rowub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. Equality constraints are specified by setting the row lower and upper bounds to the same value.

A variant usage of this function is to consider the linear programming problem in the following form,

    linprog(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

```math
\\begin{align*}
\\text{min}_{x} \\quad & c^{T}x \\\\
\\text{s.t.}    \\quad & A_{i}^{T} \\, \\text{sense}_{i} \\, b_{i} \\quad \\forall i \\\\
                \\quad & l \\leq x \\leq u
\\end{align*}
```

where:

* ``c``: is the objective vector, always in the sense of minimization

* ``A``: is the constraint matrix, with rows a_i (viewed as column-oriented vectors)

* sense: is a vector of constraint sense characters ``<``, ``=``, and ``>``

* ``b``: is the right-hand side vector

* ``l``: is the vector of lower bounds on the variables

* ``u`` is the vector of upper bounds on the variables, and solver specifies the desired solver, see [Choosing Solvers](@ref).


A shortened version is defined as::

    linprog(c, A, lb, ub, solver) = linprog(c, A, lb, ub, 0, Inf, solver)

!!! note
    The function [`linprog`](@ref) calls two independent functions for building and solving the linear programming problem, namely [`buildlp`](@ref) and [`solvelp`](@ref).

The [`linprog`](@ref) function returns an instance of the type::

```julia
    type LinprogSolution
        status
        objval
        sol
        attrs
    end
```

where `status` is a termination status symbol, one of `:Optimal`, `:Infeasible`, `:Unbounded`, `:UserLimit` (iteration limit or timeout), `:Error` (and maybe others).

If `status` is `:Optimal`, the other members have the following values:

* `objval` -- optimal objective value
* `sol` -- primal solution vector
* `attrs` -- a dictionary that may contain other relevant attributes such as:

  - `redcost` -- dual multipliers for active variable bounds (zero if inactive)
  - `lambda` -- dual multipliers for active linear constraints (equalities are always active)

If `status` is `:Infeasible`, the `attrs` member will contain an `infeasibilityray` if available; similarly for `:Unbounded` problems, `attrs` will contain an `unboundedray` if available.

  - `colbasis` -- optimal simplex basis statuses for the variables (columns) if available. Possible values are `:NonbasicAtLower`, `:NonbasicAtUpper`, `:Basic`, and `:Superbasic` (not yet implemented by any solvers)
  - `rowbasis` -- optimal simplex basis statuses for the constraints (rows) if available (not yet implemented by any solvers)

For example, we can solve the two-dimensional problem (see ``test/linprog.jl``):

```math
    \\begin{align*}
    \\text{min}_{x,y} \\quad & -x \\\\
    \\text{s.t.}      \\quad & 2x + y \\leq 1.5 \\\\
                      \\quad & x \\geq 0, y \\geq 0
    \\end{align*}
```

```julia
    using MathProgBase, Clp

    sol = linprog([-1,0],[2 1],'<',1.5, ClpSolver())
    if sol.status == :Optimal
        println("Optimal objective value is \$(sol.objval)")
        println("Optimal solution vector is: [\$(sol.sol[1]), \$(sol.sol[2])]")
    else
        println("Error: solution status \$(sol.status)")
    end
```

"""
function linprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)
    m = buildlp(c, A, rowlb, rowub, lb, ub, solver)
    return solvelp(m)
end

linprog(c,A,rowlb,rowub, solver::AbstractMathProgSolver) = linprog(c,A,rowlb,rowub,0,Inf, solver)

buildlp(c,A,rowlb,rowub, solver::AbstractMathProgSolver) = buildlp(c,A,rowlb,rowub,0,Inf, solver)

# Old versions
linprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector) = no_lp_solver()
linprog(c,A,rowlb,rowub) = no_lp_solver()

export linprog, buildlp, solvelp
