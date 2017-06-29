
type QuadprogSolution
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

"""

    quadprog(c::InputVector, Q::AbstractMatrix, A::AbstractMatrix, sense::InputVector, b::InputVector, lb::InputVector, ub::InputVector, solver::AbstractMathProgSolver)

Solves the quadratic programming problem:

```math
    \\begin{align*}
    \\text{min}_{x} \\quad & \\frac{1}{2}x^TQx + c^Tx \\\\
    \\text{s.t.}    \\quad & a_i^Tx \\, \\text{sense}_i \\, b_i \\forall \\quad i \\\\
                    \\quad & l \\leq x \\leq u \\\\
    \\end{align*}
```

where:

*    ``c`` is the objective vector, always in the sense of minimization

*    ``Q`` is the Hessian matrix of the objective

*    ``A`` is the constraint matrix, with rows :math:`a_i` (viewed as column-oriented vectors)

*    ``sense`` is a vector of constraint sense characters ``<``, ``=``, and ``>``

*    ``b`` is the right-hand side vector

*    ``l`` is the vector of lower bounds on the variables

*    ``u`` is the vector of upper bounds on the variables, and

*    ``solver`` specifies the desired solver, see [Choosing Solvers](@ref).

A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. note::
    Quadratic programming solvers extensively exploit the sparsity of the Hessian matrix ``Q`` and the constraint matrix ``A``. While both dense and sparse matrices are accepted, for large-scale problems sparse matrices should be provided if permitted by the problem structure.

The [`quadprog`](@ref) function returns an instance of the type::

```julia
type QuadprogSolution
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
* `attrs` -- a dictionary that may contain other relevant attributes (not currently used).

Analogous shortened and range-constraint versions are available as well.

We can solve the three-dimensional QP (see `test/quadprog.jl`):

```math
    \\begin{align*}
    \\text{min}_{x,y,z} \\quad & x^2+y^2+z^2+xy+yz \\\\
    \\text{s.t.}        \\quad & x + 2y + 3z \\geq 4 \\\\
                        \\quad & x + y \\geq 1
    \\end{align*}
```

```julia
using MathProgBase, Ipopt

sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf, IpoptSolver())
if sol.status == :Optimal
    println("Optimal objective value is \$(sol.objval)")
    println("Optimal solution vector is: [\$(sol.sol[1]), \$(sol.sol[2]), \$(sol.sol[3])]")
else
    println("Error: solution status \$(sol.status)")
end
```

"""
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
