abstract MathProgCallbackData
export MathProgCallbackData

# callback has signature:
# function callback(d::MathProgCallbackData)

@define_interface begin
    setlazycallback!
    setcutcallback!
    setheuristiccallback!
    cbgetmipsolution
    cbgetlpsolution
    cbgetobj
    cbgetbestbound
    cbgetexplorednodes
    cbgetstate
    cbaddcut!
    cbaddlazy!
    cbaddsolution!
    cbsetsolutionvalue!
end

