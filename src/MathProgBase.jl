VERSION >= v"0.4.0-dev+6521" && __precompile__()

module MathProgBase

include(joinpath(dirname(@__FILE__),"SolverInterface","SolverInterface.jl"))
using .SolverInterface

# deprecated name
if VERSION >= v"0.4-rc2"
    @eval @Base.deprecate_binding MathProgSolverInterface SolverInterface
else
    const MathProgSolverInterface = SolverInterface
end

include("defaultsolvers.jl")

include(joinpath(dirname(@__FILE__),"HighLevelInterface","HighLevelInterface.jl"))
using .HighLevelInterface
export linprog, mixintprog, quadprog

end
