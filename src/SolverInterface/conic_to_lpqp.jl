# wrapper to convert Conic solver into LPQP solver

# To enable LPQP support from a Conic solver, define, e.g.,
# LinearQuadraticModel(s::ECOSSolver) = ConicToLPQPBridge(ConicModel(s))

mutable struct ConicToLPQPBridge <: AbstractLinearQuadraticModel
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
    RSOCconstrs::Vector{Vector{Int}} # x'x <= y*z, y is first index in vector, z is second
    vartypes::Vector{Symbol}
end

ConicToLPQPBridge(s::AbstractConicModel) = ConicToLPQPBridge(s, sparse(Int[],Int[],Float64[]), Float64[], Float64[], Float64[], Float64[], Float64[], :Min, Int[], Array{Vector{Int}}(undef, 0),Array{Vector{Int}}(undef, 0), Symbol[])

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

const notsoc_error = "For conic solvers, only quadratic constraints in second-order cone format (x'x <= y^2) or rotated second-order cone format (x'x <= yz) are supported"

function addquadconstr!(wrap::ConicToLPQPBridge, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    if length(linearidx) > 0 || length(linearval) > 0 || sense != '<' || rhs != 0 || count((quadval .== 0.0) .| (quadval .== 1.0)) != length(quadval) - 1 || count(quadval .== -1.0) != 1
        error(notsoc_error)
    end
    length(quadrowidx) == length(quadcolidx) == length(quadval) || error("Inconsistent dimensions")
    # check if SOC or RSOC or neither
    n_pos_on_diag = 0
    off_diag_idx  = 0
    neg_diag_idx  = 0
    n = length(quadval)
    nz = 0
    for i in 1:n
        quadval[i] == 0 && continue
        nz += 1
        if quadrowidx[i] == quadcolidx[i]
            if quadval[i] == 1
                n_pos_on_diag += 1
            elseif quadval[i] == -1
                neg_diag_idx = i
            else
                error(notsoc_error)
            end
        else
            if quadval[i] == -1
                if !(neg_diag_idx == off_diag_idx == 0)
                    error(notsoc_error)
                end
                off_diag_idx = i
            else
                error(notsoc_error)
            end
        end
    end
    if n_pos_on_diag == nz-1 && neg_diag_idx > 0
        # SOC
        SOCconstr = Int[]
        push!(SOCconstr, quadcolidx[neg_diag_idx])
        for i in 1:n
            if quadval[i] == 1
                push!(SOCconstr,quadcolidx[i])
            end
        end
        push!(wrap.SOCconstrs, SOCconstr)
    elseif n_pos_on_diag == nz-1 && off_diag_idx > 0
        # Rotated SOC
        RSOCconstr = Int[]
        push!(RSOCconstr, quadcolidx[off_diag_idx])
        push!(RSOCconstr, quadrowidx[off_diag_idx])
        for i in 1:n
            if quadval[i] == 1
                push!(RSOCconstr,quadcolidx[i])
            end
        end
        push!(wrap.RSOCconstrs, RSOCconstr)

    else
        error(notsoc_error)
    end
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

    constr_cones = []
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
            push!(constr_cones, (:Zero,extrarows:extrarows))
            push!(wrap.varboundmap, j)
            continue
        end
        if collb[j] != -Inf
            # Variable has lower bound
            extrarows += 1
            push!(I,extrarows)
            push!(J,j)
            push!(b, collb[j])
            push!(constr_cones, (:NonPos,extrarows:extrarows))
            push!(wrap.varboundmap, j)
        end
        if colub[j] != Inf
            # Variable has upper bound
            extrarows += 1
            push!(I,extrarows)
            push!(J,j)
            push!(b, colub[j])
            push!(constr_cones, (:NonNeg,extrarows:extrarows))
            push!(wrap.varboundmap, j)
        end
    end
    #@show extrarows

    for it in 1:nrow
        # Equality constraint
        if rowlb[it] == rowub[it]
            # a'x = b ==> b - a'x = 0
            push!(b, rowlb[it])
            push!(constr_cones,(:Zero,it.+(extrarows:extrarows)))
        # Range constraint - not supported
        elseif rowlb[it] != -Inf && rowub[it] != Inf
            error("Ranged constraints unsupported!")
        # Less-than constraint
        elseif rowlb[it] == -Inf
            # a'x <= b ==> b - a'x >= 0
            push!(b, rowub[it])
            push!(constr_cones,(:NonNeg,it.+(extrarows:extrarows)))
        # Greater-than constraint
        else
            # a'x >= b ==> b - a'x <= 0
            push!(b, rowlb[it])
            push!(constr_cones,(:NonPos,it.+(extrarows:extrarows)))
        end
    end

    # Now SOC constraints.
    # We can't trivially apply these to variables because they might have overlapping indices.
    # Instead, append to A.
    ISOC = Int[]
    JSOC = Int[]
    VSOC = Float64[]
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
            push!(VSOC, -1.0)
            push!(b, 0.0)
        end
        push!(constr_cones, (:SOC, (nrow+extrarows+soc_row+1):(nrow+extrarows+soc_row+n)))
        soc_row += n
    end
    RSOCconstrs = wrap.RSOCconstrs
    for it in 1:length(RSOCconstrs)
        n = length(RSOCconstrs[it])
        if collb[RSOCconstrs[it][1]] < 0 || collb[RSOCconstrs[it][2]] < 0
            error("Invalid rotated second-order conic constraint: x'x <= yz requires y >= 0 and z >= 0")
        end
        for k in 1:n
            push!(ISOC, soc_row+k)
            push!(JSOC, RSOCconstrs[it][k])
            if k <= 2
                push!(VSOC, -1.0/sqrt(2))
            else
                push!(VSOC, -1.0)
            end
            push!(b, 0.0)
        end
        push!(constr_cones, (:SOCRotated, (nrow+extrarows+soc_row+1):(nrow+extrarows+soc_row+n)))
        soc_row += n
    end


    A = vcat(sparse(I,J,ones(extrarows),extrarows, nvar),
             A,
             sparse(ISOC,JSOC,VSOC,soc_row,nvar))

    if wrap.sense == :Max
        obj = -obj
    end

    #@show obj, full(A), b, constr_cones, var_cones
    loadproblem!(wrap.m, obj, A, b, constr_cones, var_cones)
    if !all(t -> t == :Cont, wrap.vartypes)
        setvartype!(wrap.m, wrap.vartypes)
    end
    optimize!(wrap.m)

end

getsolution(wrap::ConicToLPQPBridge) = getsolution(wrap.m)
function getconstrsolution(wrap::ConicToLPQPBridge)
    wrap.A * getsolution(wrap.m)
end
status(wrap::ConicToLPQPBridge) = status(wrap.m)
for f in [:getobjval, :getobjbound]
    @eval function ($f)(wrap::ConicToLPQPBridge)
        if wrap.sense == :Max
            return -$f(wrap.m)
        else
            return $f(wrap.m)
        end
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
        Compat.rmul!(redcost, -1.0)
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
        Compat.rmul!(constrduals, -1.0)
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
function addvar!(wrap::ConicToLPQPBridge, constridx::AbstractArray{T}, constrcoef, l, u, objcoef) where T<:Integer
    wrap.A = [wrap.A sparsevec(constridx, constrcoef, size(wrap.A, 1))]
    push!(wrap.collb, l)
    push!(wrap.colub, u)
    push!(wrap.obj, objcoef)
end
function addconstr!(wrap::ConicToLPQPBridge, varidx::AbstractArray{T}, coef, lb, ub) where T<:Integer
    wrap.A = [wrap.A; sparsevec(varidx, coef, size(wrap.A, 2))']
    push!(wrap.rowlb, lb)
    push!(wrap.rowub, ub)
end
function setvartype!(wrap::ConicToLPQPBridge, v)
    wrap.vartypes = copy(v)
end
