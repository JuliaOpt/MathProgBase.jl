
mutable struct MixintprogSolution
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
        rowlb = Array{realtype}(undef, nrow)
        rowub = Array{realtype}(undef, nrow)
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


