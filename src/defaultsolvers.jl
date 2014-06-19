# these are the default solver packages in order of decreasing precedence
# for each supported problem class.

const LPsolvers = [(:Clp,:ClpSolver),
                   (:GLPKMathProgInterface,:GLPKSolverLP),
                   (:Gurobi,:GurobiSolver),
                   (:CPLEX,:CplexSolver),
                   (:Mosek,:MosekSolver)]

const MIPsolvers = [(:Cbc,:CbcSolver),
                    (:GLPKMathProgInterface,:GLPKSolverMIP),
                    (:Gurobi,:GurobiSolver),
                    (:CPLEX,:CplexSolver),
                    (:Mosek,:MosekSolver)]
     
const QPsolvers = [(:Gurobi,:GurobiSolver),
                   (:CPLEX,:CplexSolver),
                   (:Mosek,:MosekSolver)]

const SDPsolvers = [(:Mosek,:MosekSolver)]

using Base.Meta

# Don't load packages for default solvers until needed.
# This reduces the startup time for MathProgBase.
for solvertype in ["LP", "MIP", "QP", "SDP"]
    typename = symbol("Default"*solvertype*"Solver")
    @eval begin
        type $typename <: MathProgSolverInterface.AbstractMathProgSolver
        end
    end
    defaultname = symbol("default"*solvertype*"solver")
    @eval $defaultname = ($typename)()

    solvers = symbol(solvertype*"solvers")


    @eval begin function MathProgSolverInterface.model(s::$typename)
        for (pkgname, solvername) in $solvers
            if isdir(Pkg.dir((string(pkgname))))
                try
                    eval(Expr(:import,pkgname))
                catch
                    warn("Package ",string(pkgname)," is installed but couldn't be loaded")
                end
                ex = Expr(:(=), $(quot(defaultname)), Expr(:call,Expr(:.,pkgname,quot(solvername))))
                eval(ex)
                ex = Expr(:call,:model, $(quot(defaultname)))
                return eval(ex)
            end
        end
        suggestions = join(["\"$(pkgname)\", " for (pkgname,solvername) in $solvers], ' ')
        error("No ",$solvertype, " solver detected. Try installing one of the following packages: ", suggestions, " and restarting Julia")
    end end
end
