module MathProgBase

require(joinpath(Pkg.dir("MathProgBase"),"src","SolverInterface","SolverInterface.jl"))
using SolverInterface

# deprecated name
const MathProgSolverInterface = SolverInterface

include("defaultsolvers.jl")

require(joinpath(Pkg.dir("MathProgBase"),"src","HighLevelInterface","HighLevelInterface.jl"))
using HighLevelInterface
export linprog, mixintprog, quadprog

end
