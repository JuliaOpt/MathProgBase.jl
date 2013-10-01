module MathProgBase

    # import all available solvers 
    const solvers = [:Clp, :Cbc, :GLPKMathProgInterface, :Gurobi]
    for s in solvers
        if Pkg.installed(string(s)) != nothing
            eval(Expr(:import,s))
        end
    end

    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    using LinprogSolverInterface
    include("linprog.jl")
    include("mixintprog.jl")
end
