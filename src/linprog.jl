if Pkg.installed("Clp") != nothing
    lpsolver = Clp.ClpSolverInterface
elseif Pkg.installed("GLPKMathProgInterface") != nothing
    lpsolver = GLPKInterfaceLP
elseif Pkg.installed("Gurobi") != nothing
    lpsolver = Gurobi
else
    lpsolver = nothing
end
setlpsolver(s) = (global lpsolver; lpsolver = s)
function setlpsolver(s::Symbol)
    global lpsolver
    if s == :Clp
        lpsolver = Clp
    elseif s == :GLPK
        lpsolver = GLPKInterfaceLP
    elseif s == :Gurobi
        lpsolver = Gurobi
    else
        error("Unrecognized LP solver name $s")
    end
end

export setlpsolver

type LinprogSolution
    status
    objval
    sol
    attrs
end

typealias InputVector{T<:Real} Union(Vector{T},Real)

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



function linprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector; options...)
    if lpsolver == nothing
        error("No LP solver installed. " *
              "Please run Pkg.add(\"Clp\") or Pkg.add(\"GLPKMathProgInterface\") " *
              "and reload MathProgBase")
    end
    m = lpsolver.model(;options...)
    nrow,ncol = size(A)

    c = expandvec(c, ncol)
    rowlbtmp = expandvec(rowlb, nrow)
    rowubtmp = expandvec(rowub, nrow)
    lb = expandvec(lb, ncol)
    ub = expandvec(ub, ncol)
    
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
    optimize(m)
    stat = status(m)
    if stat == :Optimal
        attrs = Dict()
        attrs[:redcost] = getreducedcosts(m)
        attrs[:lambda] = getconstrduals(m)
        return LinprogSolution(stat, getobjval(m), getsolution(m), attrs)
    else
        return LinprogSolution(stat, nothing, [], Dict())
    end
end

linprog(c,A,rowlb,rowub; options...) = linprog(c,A,rowlb,rowub,0,Inf; options...)

export linprog


