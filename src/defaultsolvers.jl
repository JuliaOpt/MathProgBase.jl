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
                   (:Mosek,:MosekSolver),
                   (:Ipopt,:IpoptSolver)]

const SDPsolvers = [(:Mosek,:MosekSolver),
                    (:SCS,:SCSSolver)]

const NLPsolvers = [(:Ipopt,:IpoptSolver)]

const Conicsolvers = [(:ECOS,:ECOSSolver),
                      (:SCS,:SCSSolver)]

using Base.Meta

# Don't load packages for default solvers until needed.
# This reduces the startup time for MathProgBase.
for solvertype in ["LP", "MIP", "QP", "SDP", "NLP", "Conic"]
    typename = symbol("Default"*solvertype*"Solver")
    @eval begin
        type $typename <: SolverInterface.AbstractMathProgSolver
        end
    end
    defaultname = symbol("default"*solvertype*"solver")
    @eval $defaultname = ($typename)()

    solvers = symbol(solvertype*"solvers")


    for t in (:LinearQuadraticModel,:ConicModel,:NonlinearModel)
        @eval function SolverInterface.$t(s::$typename)
            for (pkgname, solvername) in $solvers
                try
                    eval(Expr(:import,pkgname))
                catch
                    if isdir(Pkg.dir((string(pkgname))))
                        warn("Package ",string(pkgname),
                             " is installed but couldn't be loaded. ",
                             "You may need to run `Pkg.build(\"$pkgname\")`")
                    end
                    continue
                end
                ex = Expr(:(=), $(quot(defaultname)),
                          Expr(:call, Expr(:., pkgname, quot(solvername))))
                eval(ex)
                ex = Expr(:call, $(quot(t)), $(quot(defaultname)))
                return eval(ex)
            end
            suggestions = join(["\"$(pkgname)\", " for (pkgname,solvername) in $solvers], ' ')
            error("No ",$solvertype, " solver detected. Try installing one of the following packages: ", suggestions, " and restarting Julia")
        end
    end
end
