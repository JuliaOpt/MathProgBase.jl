module MathProgBase

    if Pkg.installed("Clp") != nothing
        @eval import Clp
    end
    if Pkg.installed("CoinMP") != nothing
        @eval import CoinMP
    end
    if Pkg.installed("GLPKMathProgInterface") != nothing
        @eval using GLPKMathProgInterface
    end

    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    using LinprogSolverInterface
    include("linprog.jl")
    include("mixintprog.jl")
end
