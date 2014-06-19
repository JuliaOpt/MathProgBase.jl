
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
    cbaddcut!,
    cbaddlazy!,
    cbaddsolution!,
    cbsetsolutionvalue!


abstract MathProgCallbackData

# callback has signature:
# function callback(d::MathProgCallbackData)

for func in [:setlazycallback!,
             :setcutcallback!,
             :setheuristiccallback!,
             :cbgetmipsolution,
             :cbgetlpsolution,
             :cbgetobj,
             :cbgetbestbound,
             :cbgetexplorednodes,
             :cbgetstate,
             :cbaddcut!,
             :cbaddlazy!,
             :cbaddsolution!,
             :cbsetsolutionvalue!]
    @eval $(func)() = throw(MethodError($(func),()))
end
