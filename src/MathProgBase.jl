__precompile__()

module MathProgBase

include(joinpath(dirname(@__FILE__),"SolverInterface","SolverInterface.jl"))

include(joinpath(dirname(@__FILE__),"HighLevelInterface","HighLevelInterface.jl"))
export linprog, mixintprog, quadprog, buildlp, solvelp

end
