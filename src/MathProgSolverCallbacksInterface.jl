
export MathProgCallbackData,
    setlazycallback!,
    setcutcallback!,
    setheuristiccallback!,
    cbgetmipsolution,
    cbgetlpsolution,
    cbgetobj,
    cbgetbestbound,
    cbgetexplorednodes,
    cbgetstate,
    cbaddsolution!,
    cbaddcut!,
    cbaddlazy!

abstract MathProgCallbackData

# callback has signature:
# function callback(d::MathProgCallbackData)

# set to nothing to clear callback
setlazycallback!(m::AbstractMathProgModel,f) = error("Not implemented")
setcutcallback!(m::AbstractMathProgModel,f) = error("Not implemented")
setheuristiccallback!(m::AbstractMathProgModel,f) = error("Not implemented")

cbgetmipsolution(d::MathProgCallbackData) = error("Not Implemented")
cbgetlpsolution(d::MathProgCallbackData) = error("Not Implemented")
# best feasible objective value
cbgetobj(d::MathProgCallbackData) = error("Not Implemented")
cbgetbestbound(d::MathProgCallbackData) = error("Not implemented")
cbgetexplorednodes(d::MathProgCallbackData) = error("Not implemented")

# returns :MIPNode :MIPSol :Other
cbgetstate(d::MathProgCallbackData) = error("Not implemented")

cbaddsolution!(d::MathProgCallbackData,x) = error("Not implemented")
cbaddcut!(d::MathProgCallbackData,varidx,varcoef,sense,rhs) = error("Not implemented")
cbaddlazy!(d::MathProgCallbackData,varidx,varcoef,sense,rhs) = error("Not implemented")

