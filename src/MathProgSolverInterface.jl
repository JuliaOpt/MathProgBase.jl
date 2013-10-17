module MathProgSolverInterface

export MathProgSolver,
    loadproblem,
    writeproblem,
    getvarLB,
    setvarLB,
    getvarUB,
    setvarUB,
    getconstrLB,
    setconstrLB,
    getconstrUB,
    setconstrUB,
    getobj,
    setobj,
    addvar,
    addconstr,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    optimize,
    status,
    getobjval,
    getobjbound,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver,
    setvartype,
    getvartype

abstract MathProgSolver



loadproblem(m::MathProgSolver, filename::String) = error("Not Implemented")
loadproblem(m::MathProgSolver, A, collb, colub, obj, rowlb, rowub) = error("Not Implemented")

writeproblem(m::MathProgSolver, filename::String) = error("Not Implemented")

getvarLB(m::MathProgSolver) = error("Not Implemented")
setvarLB(m::MathProgSolver, collb) = error("Not Implemented")
getvarUB(m::MathProgSolver) = error("Not Implemented")
setvarUB(m::MathProgSolver, colub) = error("Not Implemented")
getconstrLB(m::MathProgSolver) = error("Not Implemented")
setconstrLB(m::MathProgSolver, rowlb) = error("Not Implemented")
getconstrUB(m::MathProgSolver) = error("Not Implemented")
setconstrUB(m::MathProgSolver, rowub) = error("Not Implemented")
getobj(m::MathProgSolver) = error("Not Implemented")
setobj(m::MathProgSolver, obj) = error("Not Implemented")

addvar(m::MathProgSolver, rowidx, rowcoef, collb, colub, objcoef) = error("Not Implemented")

addconstr(m::MathProgSolver, colidx, colcoef, rowlb, rowub) = error("Not Implemented")

updatemodel(m::MathProgSolver) = error("Not Implemented")

setsense(m::MathProgSolver,sense) = error("Not Implemented")
getsense(m::MathProgSolver) = error("Not Implemented")

numvar(m::MathProgSolver) = error("Not Implemented")
numconstr(m::MathProgSolver) = error("Not Implemented")

optimize(m::MathProgSolver) = error("Not Implemented")

status(m::MathProgSolver) = error("Not Implemented")

getobjval(m::MathProgSolver) = error("Not Implemented")

getobjbound(m::MathProgSolver) = error("Not Implemented")

getsolution(m::MathProgSolver) = error("Not Implemented")

getconstrsolution(m::MathProgSolver) = error("Not Implemented")

getreducedcosts(m::MathProgSolver) = error("Not Implemented")
getconstrduals(m::MathProgSolver) = error("Not Implemented")

getrawsolver(m::MathProgSolver) = error("Not Implemented")

setvartype(m::MathProgSolver, vartype) = error("Not Implemented")
getvartype(m::MathProgSolver) = error("Not Implemented")

end
