if Pkg.installed("CoinMP") != nothing
    mipsolver = CoinMP
elseif Pkg.installed("GLPKMathProgInterface") != nothing
    mipsolver = GLPKInterfaceMIP
elseif Pkg.installed("Gurobi") != nothing
    mipsolver = Gurobi
else
    mipsolver = nothing
end
setmipsolver(s) = (global mipsolver; mipsolver = s)
function setmipsolver(s::Symbol)
    global mipsolver
    if s == :CoinMP
        mipsolver = CoinMP
    elseif s == :GLPK
        mipsolver = GLPKInterfaceMIP
    elseif s == :Gurobi
        mipsolver = Gurobi
    else
        error("Unrecognized MIP solver name $s")
    end
end

export setmipsolver

type MixintprogSolution
    status
    objval
    sol
    attrs
end

typealias CharInputVector Union(Vector{Char},Real)

function mixintprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, vartypes::CharInputVector, lb::InputVector, ub::InputVector; options...)
    if mipsolver == nothing
        error("No MIP solver installed. " *
              "Please run Pkg.add(\"CoinMP\") or Pkg.add(\"GLPKMathProgInterface\") " *
              "and reload MathProgBase")
    end
    m = mipsolver.model(;options...)
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
    
    loadproblem(m, A, lb, ub, c, rowlb, rowub)
    setvartype(m, vartypes)
    optimize(m)
    stat = status(m)
    if stat == :Optimal
        attrs = Dict()
        attrs[:objbound] = getobjbound(m)
        return MixintprogSolution(stat, getobjval(m), getsolution(m), attrs)
    else
        return MixintprogSolution(stat, nothing, [], Dict())
    end
end

mixintprog(c,A,rowlb,rowub,varypes) = mixintprog(c,A,rowlb,rowub,vartypes,0,Inf)

export mixintprog


