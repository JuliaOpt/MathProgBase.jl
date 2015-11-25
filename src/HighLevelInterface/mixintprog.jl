
type MixintprogSolution
    status
    objval
    sol
    attrs
end

typealias SymbolInputVector Union{Vector{Symbol},Symbol}

function mixintprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, vartypes::SymbolInputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver = MathProgBase.defaultMIPsolver)
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
        rowlb = Array(realtype, nrow)
        rowub = Array(realtype, nrow)
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

mixintprog(c,A,rowlb,rowub,vartypes,solver::AbstractMathProgSolver=defaultMIPsolver) = mixintprog(c,A,rowlb,rowub,vartypes,0,Inf,solver)

export mixintprog


