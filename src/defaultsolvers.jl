# these are the default solver packages in order of decreasing precedence
# for each supported problem class.

const LPsolvers = [(:Clp,:ClpSolver),
                   (:GLPKMathProgInterface,:GLPKSolverLP),
                   (:Gurobi,:GurobiSolver),
                   (:CPLEX,:CplexSolver),
                   (:Mosek,:MosekSolver),
                   (:Xpress,:XpressSolver)]

const MIPsolvers = [(:Cbc,:CbcSolver),
                    (:GLPKMathProgInterface,:GLPKSolverMIP),
                    (:Gurobi,:GurobiSolver),
                    (:CPLEX,:CplexSolver),
                    (:Mosek,:MosekSolver),
                    (:Xpress,:XpressSolver)]

const QPsolvers = [(:Gurobi,:GurobiSolver),
                   (:CPLEX,:CplexSolver),
                   (:Mosek,:MosekSolver),
                   (:Ipopt,:IpoptSolver),
                   (:Xpress,:XpressSolver)]

const SDPsolvers = [(:Mosek,:MosekSolver),
                    (:SCS,:SCSSolver)]

const NLPsolvers = [(:Ipopt,:IpoptSolver),
                    (:KNITRO,:KnitroSolver),
                    (:Mosek,:MosekSolver)]

const Conicsolvers = [(:ECOS,:ECOSSolver),
                      (:SCS,:SCSSolver),
                      (:Mosek,:MosekSolver)]

function loaddefaultsolvers()
    for solverlist in [LPsolvers,MIPsolvers,QPsolvers,SDPsolvers,NLPsolvers,Conicsolvers]
        for (pkgname,solvername) in solverlist
            try
                eval(Expr(:import,pkgname))
                # stop when we've found a working solver in this category
                continue
            end
        end
    end
end

using Base.Meta

for solvertype in ["LP", "MIP", "QP", "SDP", "NLP", "Conic"]
    solvers = Symbol(solvertype*"solvers")
    functionname = Symbol("default"*solvertype*"solver")
    if VERSION > v"0.6.0-"
    @eval function ($functionname)()
        for (pkgname,solvername) in $solvers
            if isdefined(Main,pkgname)
                return getfield(getfield(Main,pkgname),solvername)()
            end
        end
        solvers = [String(p) for (p,s) in $solvers]
        error("No ", $solvertype, " solver detected. The recognized solver packages are: ", solvers,". One of these solvers must be installed and explicitly loaded with a \"using\" statement.")
    end
    else
    @eval function ($functionname)()
        for (pkgname,solvername) in $solvers
            alreadydefined = isdefined(Main,pkgname)
            if !alreadydefined
                try
                    eval(Expr(:import,pkgname))
                    # if we got here, package works but wasn't loaded,
                    # print warning.
                    Base.warn_once(string("The default ", $solvertype, " package is installed but not loaded. In Julia 0.6 and later, an explicit \"using ", pkgname, "\" statement will be required in order for the solver to be detected and used as a default."))
                catch
                    continue
                end
            end
            return getfield(getfield(Main,pkgname),solvername)()
        end
        suggestions = join(["\"$(pkgname)\", " for (pkgname,solvername) in $solvers], ' ')
        error("No ",$solvertype, " solver detected. Try installing one of the following packages: ", suggestions, " and restarting Julia")
    end
    end
end
