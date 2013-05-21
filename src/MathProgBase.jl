module MathProgBase
    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    using LinprogSolverInterface
    include("linprog.jl")
    include("mixintprog.jl")
end
