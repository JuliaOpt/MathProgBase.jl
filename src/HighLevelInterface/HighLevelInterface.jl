module HighLevelInterface

using ..SolverInterface
using ..MathProgBase

function warn_no_inf(T)
    if !(isinf(typemin(T)) && isinf(typemax(T)))
        Base.warn_once("Element type $T does not have an infinite value. Note that this may artifically introduce ranged (two-sided) constraints. To avoid this, consider casting the problem data to Float64.")
    end
end

include("linprog.jl")
include("mixintprog.jl")
include("quadprog.jl")

end
