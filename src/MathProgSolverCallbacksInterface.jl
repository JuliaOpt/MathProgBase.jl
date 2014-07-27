abstract MathProgCallbackData
export MathProgCallbackData

# callback has signature:
# function callback(d::MathProgCallbackData)

funcs = [:setlazycallback!,
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
for func in funcs
    @eval $(func)() = throw(MethodError($(func),()))
    eval(Expr(:export, func))
end
