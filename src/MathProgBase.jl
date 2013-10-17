module MathProgBase

    require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
    using MathProgSolverInterface
    export SolverNameAndOptions

    include("defaultsolvers.jl")
    @setdefaultLPsolver
    @setdefaultMIPsolver

    include("linprog.jl")
    include("mixintprog.jl")
end
