macro packageImport(pkgname)
    return Expr(:import,pkgname)
end

macro setDefault(pkgname, solvername)
    :(const defaultLPsolver = ($(pkgname).$(solvername))())
end

macro setMissing(typ,pkgnames)
    :(MissingSolver($(string(typ)),$pkgnames))
end

# macros to generate code to set default solver

function setdefaultLPsolver()
    solvers = [(:Clp,:ClpSolver),
               (:GLPKMathProgInterface,:GLPKSolverLP),
               (:Gurobi,:GurobiSolver),
               (:CPLEX,:CplexSolver),
               (:Mosek,:MosekSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            # try 
                eval(Expr(:import,pkgname))
                eval( :(const defaultLPsolver = ($(pkgname).$(solvername))()) )
                return nothing
            # end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    const defaultLPsolver = @setMissing("LP",pkgnames)
    nothing
end

function setdefaultMIPsolver()
    solvers = [(:Clp,:ClpSolver),
               (:GLPKMathProgInterface,:GLPKSolverLP),
               (:Gurobi,:GurobiSolver),
               (:CPLEX,:CplexSolver),
               (:Mosek,:MosekSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            try 
                eval(Expr(:import,pkgname))
                eval( :(const defaultMIPsolver = ($(pkgname).$(solvername))()) )
                return nothing
            end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    const defaultLPsolver = @setMissing("MIP",pkgnames)
    nothing
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
                return nothing
            end
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    const defaultLPsolver = @setMissing("QP",pkgnames)
    nothing
end
