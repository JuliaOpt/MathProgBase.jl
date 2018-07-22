
mutable struct QuadprogSolution
    status
    objval
    sol
    attrs
end

function no_qp_solver()
    error("""
          A quadratic programming solver must be specified as the final argument.

          You can use some of the solvers listed on the solvers table of http://www.juliaopt.org/ that has a checkmark the Linear/Quadratic column. See the solver's documentation to confirm that it supports quadratic objectives.

          A (free) example is Ipopt.jl. Once Ipopt is installed and imported via "using Ipopt", you can specify IpoptSolver() as the solver.
          Solver options are specified by using keyword arguments to IpoptSolver().
          """)
end

function quadprog(c::InputVector, Q::AbstractMatrix, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)
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
        sense = rowlbtmp
        rhs = rowubtmp
        @assert realtype <: Real
        warn_no_inf(realtype)
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

quadprog(c,Q,A,rowlb,rowub,solver::AbstractMathProgSolver) = quadprog(c,Q,A,rowlb,rowub,0,Inf,solver)

# Old versions
quadprog(c,Q,A,rowlb,rowub) = no_qp_solver()
quadprog(c::InputVector, Q::AbstractMatrix, A::AbstractMatrix, rowlb::InputVector, rowub::InputVector, lb::InputVector, ub::InputVector) = no_qp_solver()



export quadprog
