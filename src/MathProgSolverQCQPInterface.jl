
export setquadobj!,
    addquadconstr!
    

setquadobj!(m::AbstractMathProgModel,rowidx,colidx,quadval) = error("Not implemented")

addquadconstr!(m::AbstractMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs) = error("Not implemented")
