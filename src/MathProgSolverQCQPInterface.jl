
export setquadobj!,
       addquadconstr!
    
setquadobj!(m::AbstractMathProgModel,rowidx,colidx,quadval) = error("Not implemented")
setquadobj!(m::AbstractMathProgModel,Q::Matrix{Float64}) = setquadobj!(m,sparse(Q))
function setquadobj!(m::AbstractMathProgModel,Q::SparseMatrixCSC{Float64})
  if issym(Q) || istril(Q)
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
