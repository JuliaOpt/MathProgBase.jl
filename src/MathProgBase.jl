module MathProgBase

    require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
    using MathProgSolverInterface

    include("defaultsolvers.jl")

    include("linprog.jl")
    include("mixintprog.jl")
    include("quadprog.jl")
end
