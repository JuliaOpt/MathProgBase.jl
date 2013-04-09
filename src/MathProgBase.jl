module MathProgBase
    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    include("linprog.jl")
end
