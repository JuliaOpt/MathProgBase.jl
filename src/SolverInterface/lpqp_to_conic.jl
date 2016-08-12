# wrapper to convert LPQP solver into Conic solver

# To enable Conic support from an LPQP solver, define, e.g.,
# ConicModel(s::GurobiSolver) = LPQPtoConicBridge(LinearQuadraticModel(s))
# Must also implement supportedcones(). List SOC as a supported cone if it can
# be passed as a quadratic constraint.

type LPQPtoConicBridge <: AbstractConicModel
    lpqpmodel::AbstractLinearQuadraticModel
    c
    A
    b
    constr_cones
    var_cones
end

LPQPtoConicBridge(m::AbstractLinearQuadraticModel) = LPQPtoConicBridge(m, nothing, nothing, nothing, nothing, nothing)

export LPQPtoConicBridge

numvar(m::LPQPtoConicBridge) = size(m.A,2)
numconstr(m::LPQPtoConicBridge) = size(m.A,1)

# To transform Conic problems into LinearQuadratic problems
function loadproblem!(m::LPQPtoConicBridge, c, A, b, constr_cones, var_cones)
    m.c = c
    m.A = A
    m.b = b
    m.constr_cones = constr_cones
    m.var_cones = var_cones

    # Conic form        LP form
    # min  c'x          min      c'x
    #  st b-Ax ∈ K_1     st lb <= Ax <= b
    #        x ∈ K_2         l <=  x <= u

    # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated], give up.
    bad_cones = [:SDP, :ExpPrimal, :ExpDual]
    # for each SOC and SOCRotated affine constraint, we need to add auxiliary variables
    num_orig = length(c)
    num_aux = 0
    linconstr_coneidx = Int[]
    linconstr_idx = Int[]
    socconstr_idx = Int[]
    rsocconstr_idx = Int[]
    rsoc_start_idx = Int[]
    for i in 1:length(constr_cones)
        cone, idxs = constr_cones[i]
        cone in bad_cones && error("Cone type $(cone) not supported")
        if cone == :SOC
            num_aux += length(idxs)
            append!(socconstr_idx, idxs)
        elseif cone == :SOCRotated
            num_aux += length(idxs)
            push!(rsoc_start_idx, length(rsocconstr_idx)+1)
            push!(rsoc_start_idx, length(rsocconstr_idx)+2)
            append!(rsocconstr_idx, idxs)
        else
            push!(linconstr_coneidx, i)
            if isa(idxs,Number)
                push!(linconstr_idx, idxs)
            else
                append!(linconstr_idx, idxs)
            end
        end
    end
    @assert num_aux == length(socconstr_idx) + length(rsocconstr_idx)
    for (cone,idxs) in var_cones
        cone in bad_cones && error("Cone type $(cone) not supported")
    end

    c = vcat(c,zeros(num_aux))
    # Variable bounds
    l = fill(-Inf, length(c))
    u = fill(Inf, length(c))
    for (cone,idxs) in var_cones
        if cone == :SOC
            l[idxs[1]] = 0
        elseif cone == :SOCRotated
            l[idxs[1]] = 0
            l[idxs[2]] = 0
        else
            cone_l = (cone == :Free || cone == :NonPos) ? -Inf : 0.0
            cone_u = (cone == :Free || cone == :NonNeg) ?  Inf : 0.0
            for idx in idxs
                l[idx] = cone_l
                u[idx] = cone_u
            end
        end
    end

    # set bounds for auxiliary variables
    ksoc = 1
    krsoc = 1
    for (cone,idx) in constr_cones
        if cone == :SOC
            l[num_orig + ksoc] = 0
            ksoc += length(idx)
        elseif cone == :SOCRotated
            l[num_orig + length(socconstr_idx) + krsoc] = 0
            l[num_orig + length(socconstr_idx) + krsoc + 1] = 0
            krsoc += length(idx)
        end
    end

    # matrix for linear constraints
    Alin = A[linconstr_idx,:]


    # Linear constraint bounds
    lb = Array(Float64,length(linconstr_idx))
    ub = Array(Float64,length(linconstr_idx))
    k = 1
    for (cone,idxs) in constr_cones
        if cone != :SOC && cone != :SOCRotated
            # :Zero         b - Ax = s == 0 -> Ax == b
            # :NonPos       b - Ax = s <= 0 -> Ax >= b
            # :NonNeg       b - Ax = s >= 0 -> Ax <= b
            # :Free         b - Ax = s free ->  free
            for idx in idxs
                lb[k] = (cone == :Zero || cone == :NonPos) ? b[idx] : -Inf
                ub[k] = (cone == :Zero || cone == :NonNeg) ? b[idx] :  Inf
                k += 1
            end
        end
    end

    if length(socconstr_idx) > 0
        Aaux = A[socconstr_idx,:]
        # linear constraints for aux variables
        # for each ||b - Ax|| <= c - d^Tx,
        # introduce y = b - Ax, z = c-d^Tx, and say
        # y^Ty <= z^2.
        # Ax + y = b, so we just need to append some identity columns
        Alin = [ Alin spzeros(length(linconstr_idx),length(socconstr_idx))
        Aaux speye(length(socconstr_idx)) ]
        lbaux = b[socconstr_idx]
        ubaux = lbaux
        lb = [lb; lbaux]
        ub = [ub; ubaux]
    end

    if length(rsocconstr_idx) > 0
        Aaux = A[rsocconstr_idx,:]
        # same thing for rotated SOC
        # note we adjust for factor of 2:
        # rsoc has x'x <= 2pq
        # quadratic form is x'x <= pq
        diagvec = ones(length(rsocconstr_idx))
        diagvec[rsoc_start_idx] = 1/sqrt(2)
        Alin = [ Alin spzeros(size(Alin,1),length(rsocconstr_idx))
        [Aaux spzeros(size(Aaux,1),length(socconstr_idx))] spdiagm(diagvec) ]
        lbaux = b[rsocconstr_idx]
        ubaux = lbaux
        lb = [lb; lbaux]
        ub = [ub; ubaux]
    end

    loadproblem!(m.lpqpmodel, Alin, l, u, c, lb, ub, :Min)

    # Add conic constraints

    for (cone, idx) in var_cones
        if cone == :SOC
            addquadconstr!(m.lpqpmodel, Int[], Float64[], vcat(idx), vcat(idx), [-1.0; ones(length(idx)-1)], '<', 0.0)
        elseif cone == :SOCRotated
            idx1 = idx[1]
            idx2 = idx[2]
            idxrest = idx[3:end]
            addquadconstr!(m.lpqpmodel, Int[], Float64[], [idx1;idxrest], [idx2;idxrest], [-2.0; ones(length(idxrest))], '<', 0.0)
        end
    end

    ksoc = 1
    krsoc = 1
    for (cone,idx) in constr_cones
        if cone == :SOC
            idx1 = num_orig + ksoc
            idxrest = (num_orig+ksoc+1):(num_orig+ksoc+length(idx)-1)
            addquadconstr!(m.lpqpmodel, Int[], Float64[], [idx1; idxrest], [idx1; idxrest], [-1.0;ones(length(idxrest))], '<', 0.0)
            ksoc += length(idx)
        elseif cone == :SOCRotated
            idx1 = num_orig + length(socconstr_idx) + krsoc
            idx2 = num_orig + length(socconstr_idx) + krsoc + 1
            idxrest = (idx2+1):(idx1+length(idx)-1)
            addquadconstr!(m.lpqpmodel, Int[], Float64[], [idx1; idxrest], [idx2; idxrest], [-1.0;ones(length(idxrest))], '<', 0.0)
            krsoc += length(idx)
        end
    end
end

for f in [:optimize!, :status, :getsolution, :getobjval, :getvartype]
    @eval $f(model::LPQPtoConicBridge) = $f(model.lpqpmodel)
end

setvartype!(model::LPQPtoConicBridge, vtype) = setvartype!(model.lpqpmodel, vtype)

for f in methods_by_tag[:rewrap]
    @eval $f(model::LPQPtoConicBridge) = $f(model.lpqpmodel)
end
