export addSDPvar!,
       addSDPmatrix!,
       addSDPconstr!,
       setSDPobj!,
       getSDPvarsolution

addSDPvar!(m::AbstractMathProgModel, dim) = error("Not implemented")
addSDPmatrix!(m::AbstractMathProgModel, mat) = error("Not implemented")

addSDPconstr!(m::AbstractMathProgModel, matvaridx, matcoefidx, scalidx, scalcoef, lb, ub) = error("Not implemented")

setSDPobj!(m::AbstractMathProgModel, matvaridx, matcoefidx) = error("Not implemented")

getSDPvarsolution(m::AbstractMathProgModel, idx) = error("Not implemented")
