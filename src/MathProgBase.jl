module MathProgBase

    if Pkg.installed("Clp") != nothing
        @eval import Clp
    end
    if Pkg.installed("Cbc") != nothing
        @eval import Cbc
    end
    if Pkg.installed("GLPKMathProgInterface") != nothing
        @eval using GLPKMathProgInterface
    end
    if Pkg.installed("Gurobi") != nothing
        @eval using Gurobi
    end

    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    using LinprogSolverInterface
    include("linprog.jl")
    include("mixintprog.jl")
end
