@define_interface begin
    setvarLB!
    setvarUB!
    setconstrLB!
    setconstrUB!
    setobj!
    addvar!
    addconstr!
    setsense!
    setvartype!
    setwarmstart!
    addsos1!
    addsos2!
end

# default addvar!, not adding to any existing constraints
addvar!(m::AbstractMathProgModel, collb, colub, objcoef) = addvar!(m, [], [], collb, colub, objcoef)

