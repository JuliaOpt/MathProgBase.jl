type MIPSolver
    solvermodule::Module
    options
end

function MIPSolver(s::Symbol;options...)
    local mod
    if s == :Cbc
        mod = Cbc
    elseif s == :GLPK
        mod = GLPKMathProgInterface.GLPKInterfaceMIP
    elseif s == :Gurobi
        mod = Gurobi
    else
        error("Unrecognized MIP solver name $s")
    end
    return MIPSolver(mod,options)
end

model(s::MIPSolver) = s.solvermodule.model(;s.options...)

# This is the default solver order. Free solvers first.
const defaultmipsolvers = [:Cbc, :GLPK, :Gurobi] 
# Default MIP Solver
function MIPSolver()
    for s in defaultmipsolvers
        try
            return MIPSolver(s)
        catch
            continue
        end
    end
    error("No recognized MIP solvers installed." *
          "Please run Pkg.add(\"Cbc\") or Pkg.add(\"GLPKMathProgInterface\") " *
          "and restart Julia")
end

export MIPSolver

type MixintprogSolution
    status
    objval
    sol
    attrs
end

typealias CharInputVector Union(Vector{Char},Real)

function mixintprog(c::InputVector, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, vartypes::CharInputVector, lb::InputVector, ub::InputVector, solver::MIPSolver = MIPSolver())
    m = solver.solvermodule.model(;solver.options...)
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

mixintprog(c,A,rowlb,rowub,varypes,solver::MIPSolver=MIPSolver()) = mixintprog(c,A,rowlb,rowub,vartypes,0,Inf,solver)

export mixintprog


