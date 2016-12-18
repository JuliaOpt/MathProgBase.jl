# wrapper to convert Conic solver into LPQP solver

# To enable LPQP support from a Conic solver, define, e.g.,
# LinearQuadraticModel(s::ECOSSolver) = ConicToLPQPBridge(ConicModel(s))

type ConicToLPQPBridge <: AbstractLinearQuadraticModel
    m::AbstractConicModel
    A::SparseMatrixCSC{Float64,Int}
    collb::Vector{Float64}
    colub::Vector{Float64}
    obj::Vector{Float64}
    rowlb::Vector{Float64}
    rowub::Vector{Float64}
    sense::Symbol
    varboundmap::Vector{Int}
    SOCconstrs::Vector{Vector{Int}} # x'x <= y^2, y is first index in vector
end

ConicToLPQPBridge(s::AbstractConicModel) = ConicToLPQPBridge(s, sparse(Int[],Int[],Float64[]), Float64[], Float64[], Float64[], Float64[], Float64[], :Min, Int[], Array(Vector{Int},0))

export ConicToLPQPBridge

# Loads the provided problem data to set up the linear programming problem:
# min c'x
# st  lb <= Ax <= ub
#      l <=  x <= u
# where sense = :Min or :Max
function loadproblem!(wrap::ConicToLPQPBridge, A, collb, colub, obj, rowlb, rowub, sense)
    wrap.A = A
    wrap.collb = collb
    wrap.colub = colub
    wrap.obj = obj
    wrap.rowlb = rowlb
    wrap.rowub = rowub
    wrap.sense = sense
end

const notsoc_error = "For conic solvers, only quadratic constraints in second-order conic format (x'x <= y^2) are supported"

function addquadconstr!(wrap::ConicToLPQPBridge, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    if length(linearidx) > 0 || length(linearval) > 0 || sense != '<' || rhs != 0 || mapreduce(v-> (v == 0.0 || v == 1.0), +, 0, quadval) != length(quadval) - 1 || mapreduce(v->v == -1.0, +, 0, quadval) != 1
        error(notsoc_error)
    end
    length(quadrowidx) == length(quadcolidx) == length(quadval) || error("Inconsistent dimensions")
    SOCconstr = Int[]
    for i in 1:length(quadval)
        quadrowidx[i] == quadcolidx[i] || quadval[i] == 0 || error(notsoc_error)
        if quadval[i] == 1.0
            push!(SOCconstr, quadcolidx[i])
        elseif quadval[i] == -1.0
            unshift!(SOCconstr, quadcolidx[i])
        else
            @assert quadval[i] == 0.0
        end
    end
    push!(wrap.SOCconstrs, SOCconstr)
end


function optimize!(wrap::ConicToLPQPBridge)
    A = wrap.A
    collb = wrap.collb
    colub = wrap.colub
    obj = wrap.obj
    rowlb = wrap.rowlb
    rowub = wrap.rowub
    empty!(wrap.varboundmap)
    #@show full(A), collb, colub, obj, rowlb, rowub

    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || error("Unequal lengths for row bounds")

    constr_cones = Array(Any,0)
    var_cones = [(:Free,1:nvar)]

    # for each variable bound, create a new constraint
    # x <= u  ==> u - x >= 0
    # x >= l  ==> l - x <= 0
    b = Float64[]
    I = Int[] # sparse tuples
    J = Int[]
    extrarows = 0
    for j in 1:nvar
        if collb[j] == colub[j] # Fixed variable
            extrarows += 1
            push!(I,extrarows)
            push!(J,j)
            push!(b,collb[j])
            push!(constr_cones, (:Zero,extrarows))
            push!(wrap.varboundmap, j)
            continue
        end
        if collb[j] != -Inf
            # Variable has lower bound
            extrarows += 1
            push!(I,extrarows)
            push!(J,j)
            push!(b, collb[j])
            push!(constr_cones, (:NonPos,extrarows))
            push!(wrap.varboundmap, j)
        end
        if colub[j] != Inf
            # Variable has upper bound
            extrarows += 1
            push!(I,extrarows)
            push!(J,j)
            push!(b, colub[j])
            push!(constr_cones, (:NonNeg,extrarows))
            push!(wrap.varboundmap, j)
        end
    end
    #@show extrarows

    for it in 1:nrow
        # Equality constraint
        if rowlb[it] == rowub[it]
            # a'x = b ==> b - a'x = 0
            push!(b, rowlb[it])
            push!(constr_cones,(:Zero,it+extrarows))
        # Range constraint - not supported
        elseif rowlb[it] != -Inf && rowub[it] != Inf
            error("Ranged constraints unsupported!")
        # Less-than constraint
        elseif rowlb[it] == -Inf
            # a'x <= b ==> b - a'x >= 0
            push!(b, rowub[it])
            push!(constr_cones,(:NonNeg,it+extrarows))
        # Greater-than constraint
        else
            # a'x >= b ==> b - a'x <= 0
            push!(b, rowlb[it])
            push!(constr_cones,(:NonPos,it+extrarows))
        end
    end

    # Now SOC constraints.
    # We can't trivially apply these to variables because they might have overlapping indices.
    # Instead, append to A.
    ISOC = Int[]
    JSOC = Int[]
    SOCconstrs = wrap.SOCconstrs
    soc_row = 0
    for it in 1:length(SOCconstrs)
        n = length(SOCconstrs[it])
        if collb[SOCconstrs[it][1]] < 0
            error("Invalid second-order conic constraint: x'x <= y^2 requires y >= 0")
        end
        for k in 1:n
            push!(ISOC, soc_row+k)
            push!(JSOC, SOCconstrs[it][k])
            push!(b, 0.0)
        end
        push!(constr_cones, (:SOC, (nrow+extrarows+soc_row+1):(nrow+extrarows+soc_row+n)))
        soc_row += n
    end


    A = vcat(sparse(I,J,ones(extrarows),extrarows, nvar),
             A,
             sparse(ISOC,JSOC,fill(-1.0, length(ISOC)),soc_row,nvar))

    if wrap.sense == :Max
        obj = -obj
    end

    #@show obj, full(A), b, constr_cones, var_cones
    loadproblem!(wrap.m, obj, A, b, constr_cones, var_cones)
    optimize!(wrap.m)

end

getsolution(wrap::ConicToLPQPBridge) = getsolution(wrap.m)
function getconstrsolution(wrap::ConicToLPQPBridge)
    wrap.A * getsolution(wrap.m)
end
status(wrap::ConicToLPQPBridge) = status(wrap.m)
function getobjval(wrap::ConicToLPQPBridge)
    if wrap.sense == :Max
        return -getobjval(wrap.m)
    else
        return getobjval(wrap.m)
    end
end

function getreducedcosts(wrap::ConicToLPQPBridge)
    redcost = zeros(length(wrap.collb))
    conedual = getdual(wrap.m)
    for i in 1:length(wrap.varboundmap)
        varidx = wrap.varboundmap[i]
        redcost[varidx] -= conedual[i]
    end
    if wrap.sense == :Max
        scale!(redcost,-1.0)
    end
    return redcost
end

function getconstrduals(wrap::ConicToLPQPBridge)
    constrduals = zeros(length(wrap.rowlb))
    offset = length(wrap.varboundmap)
    conedual = getdual(wrap.m)
    for i in (offset+1):(offset+length(wrap.rowlb))
        constrduals[i-offset] -= conedual[i]
    end
    if wrap.sense == :Max
        scale!(constrduals,-1.0)
    end
    return constrduals
end

numconstr(wrap::ConicToLPQPBridge) = size(wrap.A, 1)
numvar(wrap::ConicToLPQPBridge) = size(wrap.A, 2)
getvarLB(wrap::ConicToLPQPBridge) = wrap.collb
getvarUB(wrap::ConicToLPQPBridge) = wrap.colub
getconstrLB(wrap::ConicToLPQPBridge) = wrap.rowlb
getconstrUB(wrap::ConicToLPQPBridge) = wrap.rowub
getobj(wrap::ConicToLPQPBridge) = wrap.obj
getsense(wrap::ConicToLPQPBridge) = wrap.sense
function setvarLB!(wrap::ConicToLPQPBridge, l)
    wrap.collb = l
end
function setvarUB!(wrap::ConicToLPQPBridge, u)
    wrap.colub = u
end
function setconstrLB!(wrap::ConicToLPQPBridge, lb)
    wrap.rowlb = lb
end
function setconstrUB!(wrap::ConicToLPQPBridge, ub)
    wrap.rowub = ub
end
function setobj!(wrap::ConicToLPQPBridge, obj)
    wrap.obj = obj
end
function setsense!(wrap::ConicToLPQPBridge, sense)
    wrap.sense = sense
end
function addvar!{T<:Integer}(wrap::ConicToLPQPBridge, constridx::AbstractArray{T}, constrcoef, l, u, objcoef)
    wrap.A = [wrap.A sparsevec(constridx, constrcoef, size(wrap.A, 1))]
    push!(wrap.collb, l)
    push!(wrap.colub, u)
    push!(wrap.obj, objcoef)
end
function addconstr!{T<:Integer}(wrap::ConicToLPQPBridge, varidx::AbstractArray{T}, coef, lb, ub)
    wrap.A = [wrap.A; sparsevec(varidx, coef, size(wrap.A, 2))']
    push!(wrap.rowlb, lb)
    push!(wrap.rowub, ub)
end
