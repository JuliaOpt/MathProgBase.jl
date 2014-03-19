addSDPvar!(m::AbstractMathProgModel, dim) = error("Not implemented")
addSDPmatrix!(m::AbstractMathProgModel, mat) = error("Not implemented")

addSDPconstr!(m::AbstractMathProgModel, varidx, coefidx, rowlb, rowub) = error("Not implemented")

setSDPobj!(m::AbstractMathProgModel, coefidx) = error("Not implemented")

getSDPvarsolution(m::AbstractMathProgModel, idx) = error("Not implemented")
