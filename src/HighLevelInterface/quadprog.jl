
type QuadprogSolution
    status
    objval
    sol
    attrs
end

function quadprog(c::InputVector, Q::AbstractMatrix, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver = MathProgBase.defaultQPsolver)
    m = model(solver)
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
        warn_no_inf(realtype)
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
    setquadobj!(m, Q)
    optimize!(m)
    stat = status(m)
    if stat == :Optimal
        attrs = Dict()
        return QuadprogSolution(stat, getobjval(m), getsolution(m), attrs)
    else
        return QuadprogSolution(stat, nothing, [], Dict())
    end
end

quadprog(c,Q,A,rowlb,rowub,solver::AbstractMathProgSolver=defaultQPsolver) = quadprog(c,Q,A,rowlb,rowub,0,Inf,solver)

export quadprog
