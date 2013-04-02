module LinprogSolverInterface

export LinprogSolver,
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
    addcol,
    addrow,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    optimize,
    status,
    getobjval,
    getsolution,
    getrowsolution,
    getreducedcosts,
    getrowduals,
    getrawsolver

abstract LinprogSolver



loadproblem(m::LinprogSolver, filename::String) = error("Not Implemented")
loadproblem(m::LinprogSolver, A, collb, colub, obj, rowlb, rowub) = error("Not Implemented")

writeproblem(m::LinprogSolver, filename::String) = error("Not Implemented")

getvarLB(m::LinprogSolver) = error("Not Implemented")
setvarLB(m::LinprogSolver, collb) = error("Not Implemented")
getvarUB(m::LinprogSolver) = error("Not Implemented")
setvarUB(m::LinprogSolver, colub) = error("Not Implemented")
getconstrLB(m::LinprogSolver) = error("Not Implemented")
setconstrLB(m::LinprogSolver, rowlb) = error("Not Implemented")
getconstrUB(m::LinprogSolver) = error("Not Implemented")
setconstrUB(m::LinprogSolver, rowub) = error("Not Implemented")
getobj(m::LinprogSolver) = error("Not Implemented")
setobj(m::LinprogSolver, obj) = error("Not Implemented")

addcol(m::LinprogSolver, rowidx, rowcoef, collb, colub, objcoef) = error("Not Implemented")

addrow(m::LinprogSolver, colidx, colcoef, rowlb, rowub) = error("Not Implemented")

updatemodel(m::LinprogSolver) = error("Not Implemented")

setsense(m::LinprogSolver,sense) = error("Not Implemented")
getsense(m::LinprogSolver) = error("Not Implemented")

numvar(m::LinprogSolver) = error("Not Implemented")
numconstr(m::LinprogSolver) = error("Not Implemented")

optimize(m::LinprogSolver) = error("Not Implemented")

status(m::LinprogSolver) = error("Not Implemented")

getobjval(m::LinprogSolver) = error("Not Implemented")

getsolution(m::LinprogSolver) = error("Not Implemented")

getrowsolution(m::LinprogSolver) = error("Not Implemented")

getreducedcosts(m::LinprogSolver) = error("Not Implemented")
getrowduals(m::LinprogSolver) = error("Not Implemented")

getrawsolver(m::LinprogSolver) = error("Not Implemented")

end
