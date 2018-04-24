# Methods for the LinearQuadratic interface

abstract type AbstractLinearQuadraticModel <: AbstractMathProgModel end
export AbstractLinearQuadraticModel

#    writeproblem
#    freemodel!
#    optimize!

@define_interface begin
    LinearQuadraticModel
    setvarLB!
    setvarUB!
    setconstrLB!
    setconstrUB!
    setobj!
    changecoeffs!
    addvar!
    delvars!
    addconstr!
    delconstrs!
    addsos1!
    addsos2!
    writeproblem
    getvarLB
    getvarUB
    getconstrLB
    getconstrUB
    getobj
    getconstrmatrix
    numlinconstr
    getconstrsolution
    getreducedcosts
    getconstrduals
    getinfeasibilityray
    getunboundedray
    getsimplexiter, rewrap
    getbarrieriter, rewrap
    getnodecount, rewrap
    getbasis
end

# default addvar!, not adding to any existing constraints
addvar!(m::AbstractMathProgModel, collb, colub, objcoef) = addvar!(m, Int[], [], collb, colub, objcoef)

# Quadratic methods

@define_interface begin
    numquadconstr
    setquadobj!
    setquadobjterms!
    addquadconstr!
    getquadconstrsolution
    getquadconstrduals
    getquadinfeasibilityray
    getquadconstrRHS
    setquadconstrRHS!
end

function setquadobjterms!(m::AbstractLinearQuadraticModel, rowidx, colidx, quadval)
    (n = length(rowidx)) == length(colidx) == length(quadval) || error("Inconsistent input dimensions")
    nquadval = copy(quadval)
    for i in 1:length(rowidx)
        if rowidx[i] == colidx[i] # if on diagonal...
            nquadval[i] *= 2
        end
    end
    if applicable(setquadobj!, m, rowidx, colidx, nquadval)
        setquadobj!(m, rowidx, colidx, nquadval)
    else
        error("Solver does not support quadratic objectives")
    end
end

setquadobj!(m::AbstractLinearQuadraticModel,Q::Matrix) = setquadobj!(m,sparse(float(Q)))
function setquadobj!(m::AbstractLinearQuadraticModel,Q::SparseMatrixCSC{Float64})
    if issymmetric(Q) || istriu(Q)
        nnz_q = nnz(Q)
        qr = Array{Cint}(undef, nnz_q)
        qc = Array{Cint}(undef, nnz_q)
        qv = Array{Cdouble}(undef, nnz_q)
        k = 0
        colptr::Vector{Int} = Q.colptr
        nzval::Vector{Float64} = Q.nzval

        for i = 1:numvar(m)
            qi::Cint = convert(Cint, i)
            for j = colptr[i]:(colptr[i+1]-1)
                qj = convert(Cint, Q.rowval[j])
                if qi <= qj
                    k += 1
                    qr[k] = qi
                    qc[k] = qj
                    qv[k] = nzval[j]
                end
            end
        end
        if applicable(setquadobj!, m, qr, qc, qv)
            setquadobj!(m,qr[1:k],qc[1:k],qv[1:k])
        else
            error("Solver does not support quadratic objectives")
        end
    else
        error("Quadratic cost coefficient matrix Q is not symmetric or upper triangular")
    end
end
