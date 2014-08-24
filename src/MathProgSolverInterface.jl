warn("WARNING: The MathProgSolverInterface.jl file is DEPRECATED in favor of the SolverInterface/ directory. Do not point your packages directly to this file.")

require(joinpath(Pkg.dir("MathProgBase"),"src","SolverInterface","SolverInterface.jl"))
using SolverInterface

# deprecated name
const MathProgSolverInterface = SolverInterface

