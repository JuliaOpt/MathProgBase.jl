module MathProgSolverInterface

abstract AbstractMathProgModel
export AbstractMathProgModel

# immutable type which we dispatch solvers on 
abstract AbstractMathProgSolver
export AbstractMathProgSolver

# create dummy method to define function so that we can attach methods in other modules

funcs = [:model, :loadproblem!, :writeproblem,
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
         :addsos1!, :addsos2!,
         :getbasis]
for func in funcs
    @eval $(func)() = throw(MethodError($(func),()))
    eval(Expr(:export, func))
end

# default addvar!, not adding to any existing constraints
addvar!(m::AbstractMathProgModel, collb, colub, objcoef) = addvar!(m, [], [], collb, colub, objcoef)

include("MathProgSolverCallbacksInterface.jl")
include("MathProgSolverQCQPInterface.jl")
include("MathProgSolverSDPInterface.jl")
include("MathProgSolverNonlinearInterface.jl")
include("MathProgSolverConicInterface.jl")

end
