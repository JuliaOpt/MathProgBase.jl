__precompile__()

module MathProgBase

include(joinpath(dirname(@__FILE__),"SolverInterface","SolverInterface.jl"))
using .SolverInterface

include("defaultsolvers.jl")

include(joinpath(dirname(@__FILE__),"HighLevelInterface","HighLevelInterface.jl"))
using .HighLevelInterface
export linprog, mixintprog, quadprog, buildlp, solvelp

end
