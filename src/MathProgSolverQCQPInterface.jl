
export setquadobj!,
       setquadobjterms!,
       addquadconstr!

function setquadobjterms!(m::AbstractMathProgModel, rowidx, colidx, quadval)
    (n = length(rowidx)) == length(colidx) == length(quadval) || error("Inconsistent input dimensions")
    nquadval = copy(quadval)
    for i in 1:length(rowidx)
        if rowidx[i] == colidx[i] # if on diagonal...
            nquadval[i] *= 2
        end
    end
    setquadobj!(m, rowidx, colidx, nquadval)
end

setquadobj!(m::AbstractMathProgModel,rowidx,colidx,quadval) = error("Not implemented")
setquadobj!(m::AbstractMathProgModel,Q::Matrix) = setquadobj!(m,sparse(float(Q)))
function setquadobj!(m::AbstractMathProgModel,Q::SparseMatrixCSC{Float64})
  if issym(Q) || istriu(Q)
    nnz_q = nnz(Q)
    qr = Array(Cint, nnz_q)
    qc = Array(Cint, nnz_q)
    qv = Array(Cdouble, nnz_q)
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
    setquadobj!(m,qr[1:k],qc[1:k],qv[1:k])
  else
    error("Quadratic cost coefficient matrix Q is not symmetric or upper triangular")
  end
end

addquadconstr!(m::AbstractMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs) = error("Not implemented")
