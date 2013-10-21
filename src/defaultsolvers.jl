
# macros to generate code to set default solver

macro setdefaultLPsolver()
    solvers = [(:Clp,:ClpSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            importexpr = Expr(:import,pkgname)
            return esc(quote 
                $importexpr
                const defaultLPsolver = ($(pkgname).$(solvername))()
            end)
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    return esc(quote
        const defaultLPsolver = MissingSolver("LP",$pkgnames)
    end)
end

macro setdefaultMIPsolver()
    solvers = [(:Cbc,:CbcSolver)]
    for (pkgname, solvername) in solvers
        if Pkg.installed(string(pkgname)) != nothing
            importexpr = Expr(:import,pkgname)
            return esc(quote 
                $importexpr
                const defaultMIPsolver = ($(pkgname).$(solvername))()
            end)
        end
    end
    pkgnames = [pkgname for (pkgname, solvername) in solvers]
    return esc(quote
        const defaultMIPsolver = MissingSolver("MIP",$pkgnames)
    end)
end
