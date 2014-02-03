function setdefaultLPsolver()
    solvers = [(:Clp,:ClpSolver),
               (:GLPKMathProgInterface,:GLPKSolverLP),
               (:Gurobi,:GurobiSolver),
               (:CPLEX,:CplexSolver),
               (:Mosek,:MosekSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            try 
                eval(Expr(:import,pkgname))
                eval( :(const defaultLPsolver = ($(pkgname).$(solvername))()) )
                return
            catch
                warn("Package $pkgname is installed but couldn't be loaded")
            end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    eval(:(const defaultLPsolver = MissingSolver("LP",$pkgnames)))
end

function setdefaultMIPsolver()
    solvers = [(:Cbc,:CbcSolver),
               (:GLPKMathProgInterface,:GLPKSolverMIP),
               (:Gurobi,:GurobiSolver),
               (:CPLEX,:CplexSolver),
               (:Mosek,:MosekSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            try 
                eval(Expr(:import,pkgname))
                eval( :(const defaultMIPsolver = ($(pkgname).$(solvername))()) )
                return
            catch
                warn("Package $pkgname is installed but couldn't be loaded")
            end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    eval(:(const defaultMIPsolver = MissingSolver("MIP",$pkgnames)))
end

function setdefaultQPsolver()
    solvers = [(:Gurobi,:GurobiSolver),
               (:CPLEX,:CplexSolver),
               (:Mosek,:MosekSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            try 
                eval(Expr(:import,pkgname))
                eval( :(const defaultQPsolver = ($(pkgname).$(solvername))()) )
                return
            catch
                warn("Package $pkgname is installed but couldn't be loaded")
            end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    eval(:(const defaultQPsolver = MissingSolver("QP",$pkgnames)))
end
