
# addSDPvar!(m::AbstracMathProgModel)

# addSDPconstr!(m::AbstractMathProgModel, colidx, colcoef, rowlb, rowub) = error("Not Implemented")

# setSDPobj!(m::AbstractMathProgModel,rowidx,colidx,sdpval) = error("Not implemented")
# setSDPobj!(m::AbstractMathProgModel,C::Matrix) = setsdpobj!(m,sparse(float(C)))
# function setquadobj!(m::AbstractMathProgModel,C::SparseMatrixCSC{Float64})
#   if issym(C) || istriu(C)
#     nnz_q = nnz(C)
#     qr = Array(Cint, nnz_c)
#     qc = Array(Cint, nnz_c)
#     qv = Array(Cdouble, nnz_c)
#     k = 0
#     colptr::Vector{Int} = C.colptr
#     nzval::Vector{Float64} = C.nzval

#     for i = 1:numvar(m)
#       qi::Cint = convert(Cint, i)
#       for j = colptr[i]:(colptr[i+1]-1)
#         qj = convert(Cint, C.rowval[j])
#         if qi <= qj
#           k += 1
#           qr[k] = qi
#           qc[k] = qj
#           qv[k] = nzval[j]
#         end
#       end
#     end
#     setquadobj!(m,qr[1:k],qc[1:k],qv[1:k])
#   else
#     error("SDP cost coefficient matrix C is not symmetric or upper triangular")
#   end
# end

# getSDPvarsolution(m::AbstractMathProgModel, index)
