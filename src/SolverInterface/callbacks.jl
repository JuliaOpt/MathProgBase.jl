abstract type MathProgCallbackData end
export MathProgCallbackData

# callback has signature:
# function callback(d::MathProgCallbackData)

@define_interface begin
    setlazycallback!
    setcutcallback!
    setheuristiccallback!
    setinfocallback!
    cbgetmipsolution
    cbgetlpsolution
    cbgetobj
    cbgetbestbound
    cbgetexplorednodes
    cbgetstate
    cbaddcut!
    cbaddcutlocal!
    cbaddlazy!
    cbaddlazylocal!
    cbaddsolution!
    cbsetsolutionvalue!
end
