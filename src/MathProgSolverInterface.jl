module MathProgSolverInterface

export AbstractMathProgModel,
    AbstractMathProgSolver,
    MissingSolver,
    model,
    loadproblem!,
    writeproblem,
    getvarLB,
    setvarLB!,
    getvarUB,
    setvarUB!,
    getconstrLB,
    setconstrLB!,
    getconstrUB,
    setconstrUB!,
    getobj,
    setobj!,
    getconstrmatrix,
    addvar!,
    addconstr!,
    updatemodel!,
    setsense!,
    getsense,
    numvar,
    numconstr,
    optimize!,
    status,
    getobjval,
    getobjbound,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver,
    setvartype!,
    getvartype,
    getinfeasibilityray,
    getunboundedray,
    setwarmstart!,
    addsos1!,
    addsos2!

abstract AbstractMathProgModel

# immutable type which we dispatch solvers on 
abstract AbstractMathProgSolver

# create dummy method to define function so that we can attach methods in other modules

for func in [:model, :loadproblem!, :writeproblem,
             :getvarLB, :setvarLB!, :getvarUB, :setvarUB!,
             :getconstrLB, :setconstrLB!, :getconstrUB, :setconstrUB!, 
             :getobj, :setobj!, 
             :getconstrmatrix,
             :addvar!, :addconstr!,
             :updatemodel!, 
             :setsense!, :getsense,
             :numvar, :numconstr,
             :optimize!, :status,
             :getobjval, :getobjbound, 
             :getsolution, :getconstrsolution, 
             :getreducedcosts, :getconstrduals,
             :getrawsolver,
             :getvartype, :setvartype!, 
             :getinfeasibilityray, :getunboundedray,
             :setwarmstart!,
             :addsos1!, :addsos2!]
    @eval $(func)() = throw(MethodError($(func),()))
end

include("MathProgSolverCallbacksInterface.jl")
include("MathProgSolverQCQPInterface.jl")
include("MathProgSolverSDPInterface.jl")

end
