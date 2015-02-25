module MathProgBase

require(joinpath(dirname(@__FILE__),"SolverInterface","SolverInterface.jl"))
using SolverInterface

# deprecated name
const MathProgSolverInterface = SolverInterface

include("defaultsolvers.jl")

require(joinpath(dirname(@__FILE__),"HighLevelInterface","HighLevelInterface.jl"))
using HighLevelInterface
export linprog, mixintprog, quadprog

end
