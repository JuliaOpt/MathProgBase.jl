export addsdpvar!,
       addsdpmatrix!,
       addsdpconstr!,
       setsdpobj!,
       getsdpsolution

addsdpvar!(m::AbstractMathProgModel, dim) = error("Not implemented")
addsdpmatrix!(m::AbstractMathProgModel, mat) = error("Not implemented")

addsdpconstr!(m::AbstractMathProgModel, matvaridx, matcoefidx, scalidx, scalcoef, lb, ub) = error("Not implemented")

setsdpobj!(m::AbstractMathProgModel, matvaridx, matcoefidx) = error("Not implemented")

getsdpsolution(m::AbstractMathProgModel, idx) = error("Not implemented")
getsdpdual(m::AbstractMathProgModel) = errror("Not implemented")
