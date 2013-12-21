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
    setwarmstart!

abstract AbstractMathProgModel

# immutable type which we dispatch solvers on 
abstract AbstractMathProgSolver

model(s::AbstractMathProgSolver) = error("Not implemented")

immutable MissingSolver <: AbstractMathProgSolver
    solvertype::ASCIIString
    suggestions::Vector{Symbol}
end

function model(s::MissingSolver) 
    pkgnames = join(["\"$(sol)\", " for sol in s.suggestions], ' ')
    error("No $(s.solvertype) solver detected. Try installing one of the following packages: $pkgnames and restarting Julia")
end

loadproblem!(m::AbstractMathProgModel, filename::String) = error("Not Implemented")
loadproblem!(m::AbstractMathProgModel, A, collb, colub, obj, rowlb, rowub, sense) = error("Not Implemented")

writeproblem(m::AbstractMathProgModel, filename::String) = error("Not Implemented")

getvarLB(m::AbstractMathProgModel) = error("Not Implemented")
setvarLB!(m::AbstractMathProgModel, collb) = error("Not Implemented")
getvarUB(m::AbstractMathProgModel) = error("Not Implemented")
setvarUB!(m::AbstractMathProgModel, colub) = error("Not Implemented")
getconstrLB(m::AbstractMathProgModel) = error("Not Implemented")
setconstrLB!(m::AbstractMathProgModel, rowlb) = error("Not Implemented")
getconstrUB(m::AbstractMathProgModel) = error("Not Implemented")
setconstrUB!(m::AbstractMathProgModel, rowub) = error("Not Implemented")
getobj(m::AbstractMathProgModel) = error("Not Implemented")
setobj!(m::AbstractMathProgModel, obj) = error("Not Implemented")

addvar!(m::AbstractMathProgModel, rowidx, rowcoef, collb, colub, objcoef) = error("Not Implemented")

addvar!(m::AbstractMathProgModel, collb, colub, objcoef) = addvar!(m, [], [], collb, colub, objcoef)

addconstr!(m::AbstractMathProgModel, colidx, colcoef, rowlb, rowub) = error("Not Implemented")

updatemodel!(m::AbstractMathProgModel) = error("Not Implemented")

setsense!(m::AbstractMathProgModel,sense) = error("Not Implemented")
getsense(m::AbstractMathProgModel) = error("Not Implemented")

numvar(m::AbstractMathProgModel) = error("Not Implemented")
numconstr(m::AbstractMathProgModel) = error("Not Implemented")

optimize!(m::AbstractMathProgModel) = error("Not Implemented")

status(m::AbstractMathProgModel) = error("Not Implemented")

getobjval(m::AbstractMathProgModel) = error("Not Implemented")

getobjbound(m::AbstractMathProgModel) = error("Not Implemented")

getsolution(m::AbstractMathProgModel) = error("Not Implemented")

getconstrsolution(m::AbstractMathProgModel) = error("Not Implemented")

getreducedcosts(m::AbstractMathProgModel) = error("Not Implemented")
getconstrduals(m::AbstractMathProgModel) = error("Not Implemented")

getrawsolver(m::AbstractMathProgModel) = error("Not Implemented")

setvartype!(m::AbstractMathProgModel, vartype) = error("Not Implemented")
getvartype(m::AbstractMathProgModel) = error("Not Implemented")

setwarmstart!(m::AbstractMathProgModel, v) = error("Not Implemented")

include("MathProgSolverCallbacksInterface.jl")
include("MathProgSolverQCQPInterface.jl")

end
