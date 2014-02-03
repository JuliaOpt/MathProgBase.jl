module MathProgBase

    require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
    using MathProgSolverInterface
    export SolverNameAndOptions

    include("defaultsolvers.jl")
    setdefaultLPsolver()
    setdefaultMIPsolver()
    setdefaultQPsolver()

    include("linprog.jl")
    include("mixintprog.jl")
    include("quadprog.jl")
end
