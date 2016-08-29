# wrapper to convert Nonlinear solver into LinearQuadratic solver

# To enable LPQP support from a Nonlinear solver, define, e.g.,
# LinearQuadraticModel(s::IpoptSolver) = NonlinearToLPQPBridge(NonlinearModel(s))

type QuadConstr
    linearidx::Vector{Int}
    linearval::Vector{Float64}
    quadrowidx::Vector{Int}
    quadcolidx::Vector{Int}
    quadval::Vector{Float64}
    sense::Char
    rhs::Float64
end
getdata(q::QuadConstr) = q.linearidx, q.linearval, q.quadrowidx, q.quadcolidx, q.quadval

type NonlinearToLPQPBridge <: AbstractLinearQuadraticModel
    nlpmodel::AbstractNonlinearModel
    LPdata
    Qobj::Tuple{Vector{Int},Vector{Int},Vector{Float64}}
    Qconstr::Vector{QuadConstr}
    numvar::Int
    numconstr::Int
    warmstart::Vector{Float64}
    vartypes::Vector{Symbol}
end

function NonlinearToLPQPBridge(m::AbstractNonlinearModel)
    return NonlinearToLPQPBridge(m,nothing,(Int[],Int[],Float64[]),QuadConstr[],0,0,Float64[],Symbol[])
end

export NonlinearToLPQPBridge

function loadproblem!(model::NonlinearToLPQPBridge, A, l, u, c, lb, ub, sense)
    model.LPdata = (A,l,u,c,lb,ub,sense)
    model.numvar = size(A,2)
    model.numconstr = size(A,1)
end

function setquadobj!(model::NonlinearToLPQPBridge, rowidx, colidx, quadval)
    model.Qobj = (rowidx, colidx, quadval)
end

function addquadconstr!(model::NonlinearToLPQPBridge, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    push!(model.Qconstr, QuadConstr(linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs))
    model.numconstr += 1
end

function optimize!(model::NonlinearToLPQPBridge)
    A,l,u,c,lb,ub,sense = model.LPdata
    lb = copy(lb)
    ub = copy(ub)
    for q in model.Qconstr
        if q.sense == '<'
            push!(lb,-Inf)
            push!(ub,q.rhs)
        elseif q.sense == '>'
            push!(lb,q.rhs)
            push!(ub,Inf)
        else
            @assert q.sense == '='
            push!(lb,q.rhs)
            push!(ub,q.rhs)
        end
    end

    loadproblem!(model.nlpmodel, model.numvar, model.numconstr, l, u, lb, ub, sense, LPQPEvaluator(model))
    if any(model.vartypes .!= :Cont)
        setvartype!(model.nlpmodel, model.vartypes)
    end
    optimize!(model.nlpmodel)
end

for f in [:status, :getsolution, :getobjval, :getobjbound, :getreducedcosts]
    @eval $f(model::NonlinearToLPQPBridge) = $f(model.nlpmodel)
end

getsense(model::NonlinearToLPQPBridge) = model.LPdata[7]
numvar(model::NonlinearToLPQPBridge) = model.numvar
numconstr(model::NonlinearToLPQPBridge) = model.numconstr

numlinconstr(model::NonlinearToLPQPBridge) = size(model.LPdata[1],1)
numquadconstr(model::NonlinearToLPQPBridge) = length(model.Qconstr)

setvartype!(model::NonlinearToLPQPBridge,vartypes::Vector{Symbol}) = (model.vartypes = vartypes)

function getconstrsolution(model::NonlinearToLPQPBridge)
    x = getsolution(model)
    A = model.LPdata[1]
    return A*x
end

function getconstrduals(model::NonlinearToLPQPBridge)
    du = getconstrduals(model.nlpmodel)
    A = model.LPdata[1]
    return du[1:size(A,1)]
end

function getquadconstrduals(model::NonlinearToLPQPBridge)
    du = getconstrduals(model.nlpmodel)
    return du[(numlinconstr(model)+1):end]
end

type LPQPEvaluator <: AbstractNLPEvaluator
    A::SparseMatrixCSC{Float64,Int}
    c::Vector{Float64}
    Qi::Vector{Int} # quadratic objective terms
    Qj::Vector{Int}
    Qv::Vector{Float64}
    Qconstr::Vector{QuadConstr}
end

function LPQPEvaluator(model::NonlinearToLPQPBridge)
    A = model.LPdata[1]
    c = model.LPdata[4]
    Qi,Qj,Qv = model.Qobj
    return LPQPEvaluator(A,c,Qi,Qj,Qv,model.Qconstr)
end

features_available(d::LPQPEvaluator) = [:Grad, :Jac, :Hess]

function initialize(d::LPQPEvaluator, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in features_available(d))
            error("Unsupported feature $feat")
            # TODO: implement Jac-vec and Hess-vec products
            # for solvers that need them
        end
    end
end

function eval_f(d::LPQPEvaluator, x)
    obj = dot(d.c,x)
    Qi = d.Qi
    Qj = d.Qj
    Qv = d.Qv
    for k in 1:length(Qi)
        if Qi[k] == Qj[k]
            obj += Qv[k]*x[Qi[k]]*x[Qi[k]]/2
        else
            obj += Qv[k]*x[Qi[k]]*x[Qj[k]]
        end
    end
    return obj
end

function eval_grad_f(d::LPQPEvaluator, grad_f, x)
    Qi = d.Qi
    Qj = d.Qj
    Qv = d.Qv
    for j = 1:length(d.c)
        grad_f[j] = d.c[j]
    end
    for k in 1:length(Qi)
        if Qi[k] == Qj[k]
            grad_f[Qi[k]] += Qv[k]*x[Qj[k]]
        else
            grad_f[Qi[k]] += Qv[k]*x[Qj[k]]
            grad_f[Qj[k]] += Qv[k]*x[Qi[k]]
        end
    end
end

function eval_g(d::LPQPEvaluator, g, x)
    g_val = d.A*x
    g[1:size(d.A,1)] = g_val
    for i in 1:length(d.Qconstr)
        linearidx, linearval, quadrowidx, quadcolidx, quadval = getdata(d.Qconstr[i])
        val = 0.0
        for k in 1:length(linearidx)
            val += linearval[k]*x[linearidx[k]]
        end
        for k in 1:length(quadrowidx)
            val += quadval[k]*x[quadrowidx[k]]*x[quadcolidx[k]]
        end
        g[i+size(d.A,1)] = val
    end
end

function jac_structure(d::LPQPEvaluator)
    jac_nnz = nnz(d.A)
    for i in 1:length(d.Qconstr)
        jac_nnz += length(d.Qconstr[i].linearidx)
        jac_nnz += 2*length(d.Qconstr[i].quadrowidx)
    end

    I = Vector{Int}(jac_nnz)
    J = Vector{Int}(jac_nnz)
    idx = 1
    for col = 1:size(d.A,2)
        for pos = d.A.colptr[col]:(d.A.colptr[col+1]-1)
            I[idx] = d.A.rowval[pos]
            J[idx] = col
            idx += 1
        end
    end
    for i in 1:length(d.Qconstr)
        linearidx, linearval, quadrowidx, quadcolidx, quadval = getdata(d.Qconstr[i])
        for k in 1:length(linearidx)
            I[idx] = i+size(d.A,1)
            J[idx] = linearidx[k]
            idx += 1
        end
        for k in 1:length(quadrowidx)
            I[idx] = i+size(d.A,1)
            I[idx+1] = i+size(d.A,1)
            J[idx] = quadrowidx[k]
            J[idx+1] = quadcolidx[k]
            idx += 2
        end
    end
    return I,J
end

function eval_jac_g(d::LPQPEvaluator, J, x)
    idx = 1
    for col = 1:size(d.A,2)
        for pos = d.A.colptr[col]:(d.A.colptr[col+1]-1)
            J[idx] = d.A.nzval[pos]
            idx += 1
        end
    end
    for i in 1:length(d.Qconstr)
        linearidx, linearval, quadrowidx, quadcolidx, quadval = getdata(d.Qconstr[i])
        for k in 1:length(linearidx)
            J[idx] = linearval[k]
            idx += 1
        end
        for k in 1:length(quadrowidx)
            J[idx] = quadval[k]*x[quadcolidx[k]]
            J[idx+1] = quadval[k]*x[quadrowidx[k]]
            idx += 2
        end
    end
end

function hesslag_structure(d::LPQPEvaluator)
    hess_nnz = length(d.Qi)
    for i in 1:length(d.Qconstr)
        hess_nnz += length(d.Qconstr[i].quadrowidx)
    end

    I = Vector{Int}(hess_nnz)
    J = Vector{Int}(hess_nnz)

    for k in 1:length(d.Qi)
        I[k] = d.Qi[k]
        J[k] = d.Qj[k]
    end
    idx = length(d.Qi)+1
    for i in 1:length(d.Qconstr)
        linearidx, linearval, quadrowidx, quadcolidx, quadval = getdata(d.Qconstr[i])
        for k in 1:length(quadrowidx)
            qidx1 = quadrowidx[k]
            qidx2 = quadcolidx[k]
            if qidx2 > qidx1
                qidx1, qidx2 = qidx2, qidx1
            end
            I[idx] = qidx1
            J[idx] = qidx2
            idx += 1
        end
    end
    return I,J
end

function eval_hesslag(d::LPQPEvaluator, H, x, σ, μ)
    for k in 1:length(d.Qi)
        H[k] = σ*d.Qv[k]
    end
    idx = length(d.Qi) + 1
    for i in 1:length(d.Qconstr)
        linearidx, linearval, quadrowidx, quadcolidx, quadval = getdata(d.Qconstr[i])
        for k in 1:length(quadrowidx)
            l = μ[size(d.A,1)+i]
            if quadrowidx[k] == quadcolidx[k]
                H[idx] = l*2*quadval[k]
            else
                H[idx] = l*quadval[k]
            end
            idx += 1
        end
    end
end

isobjlinear(d::LPQPEvaluator) = (length(d.Qi) == 0)
isobjquadratic(d::LPQPEvaluator) = true
isconstrlinear(d::LPQPEvaluator,i::Integer) = (i <= size(d.A,1))
