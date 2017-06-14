
type MixintprogSolution
    status
    objval
    sol
    attrs
end

const SymbolInputVector = Union{Vector{Symbol},Symbol}

function no_mip_solver()
    error("""
          An integer programming solver must be specified as the final argument.

          You can use any solver listed on the solvers table of http://www.juliaopt.org/ that has checkmarks in both the Linear/Quadratic and Integer columns.

          A (free) example is Cbc.jl. Once Cbc is installed and imported via "using Cbc", you can specify CbcSolver() as the solver.
          Solver options are specified by using keyword arguments to CbcSolver().
          """)
end

"""

    mixintprog(c::InputVector, A::AbstractMatrix, sense::InputVector, b::InputVector, vartypes::SymbolInputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

Solves the same optimization problem as [`linprog`](@ref) above, except variables are additionally constrained to take only integer values if the corresponding entry in the `varypes` vector is the symbol `:Int`. Continuous variables are indicated by the value `:Cont`, binary variables should be specified by `:Bin`, semicontinuous by `:SemiCont`, and semi-integer by `:SemiInt`.

A scalar is accepted for the ``sense``, ``b``, ``vartypes``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

The [`mixintprog`](@ref) function returns an instance of the type::

    type MixintprogSolution
        status
        objval
        sol
        attrs
    end

where `status` takes the same values as with [`linprog`](@ref).

If `status` does not indicate error or infeasiblity, the other members have the following values:

* `objval` -- optimal objective value
* `sol` -- primal solution vector
* `attrs` -- a dictionary that may contain other relevant attributes such as:

  - `objbound` -- Best known lower bound on the objective value

Analogous shortened and range-constraint versions are available as well.

We can solve a [binary knapsack problem](http://en.wikipedia.org/wiki/Knapsack_problem)

```math
\\begin{align*}
    \\text{max} \\, & 5x_1 + 3x_2 + 2x_3 + 7x_4 + 4x_5 \\\\
    \\text{s.t.}    & 2x_1 + 8x_2 + 4x_3 + 2x_4 + 5x_5 \\leq 10 \\\\
                    & (x_1, x_2, x_3, x_4, x_5) \\in \\{0,1\\}^5
\\end{align*}
```

with the following code

```julia
    mixintprog(-[5.,3.,2.,7.,4.],[2. 8. 4. 2. 5.],'<',10,:Int,0,1,CbcSolver())
```

"""
function mixintprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, vartypes::SymbolInputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)
    m = LinearQuadraticModel(solver)
    nrow,ncol = size(A)

    c = expandvec(c, ncol)
    rowlbtmp = expandvec(rowlb, nrow)
    rowubtmp = expandvec(rowub, nrow)
    lb = expandvec(lb, ncol)
    ub = expandvec(ub, ncol)
    vartypes = expandvec(vartypes, ncol)

    # rowlb is allowed to be vector of senses
    if eltype(rowlbtmp) == Char
        realtype = eltype(rowubtmp)
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
    setvartype!(m, vartypes)
    optimize!(m)
    stat = status(m)
    if stat == :Optimal
        attrs = Dict()
        attrs[:objbound] = getobjbound(m)
        return MixintprogSolution(stat, getobjval(m), getsolution(m), attrs)
    else
        return MixintprogSolution(stat, nothing, [], Dict())
    end
end

mixintprog(c,A,rowlb,rowub,vartypes,solver::AbstractMathProgSolver) = mixintprog(c,A,rowlb,rowub,vartypes,0,Inf,solver)

# Old versions
mixintprog(c,A,rowlb,rowub,vartypes) = no_mip_solver()
mixintprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, vartypes::SymbolInputVector, lb::InputVector, ub::InputVector) = no_mip_solver()

export mixintprog
