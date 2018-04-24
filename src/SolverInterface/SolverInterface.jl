module SolverInterface

using Base.Meta
using Compat
using Compat.SparseArrays
using Compat.LinearAlgebra

const methods_by_tag = Dict{Symbol,Vector{Symbol}}()

macro define_interface(args)
    @assert args.head == :block
    code = quote end
    for line in args.args
        (isexpr(line, :line) || isa(line, LineNumberNode)) && continue
        (isa(line,Symbol) || isexpr(line,:tuple)) || error("Unexpected code in block")
        if isa(line,Symbol)
            fname = line
            tags = []
        else
            @assert isexpr(line,:tuple)
            fname = line.args[1]
            tags = line.args[2:end]
        end
        func = esc(fname)
        code = quote
            $code
            $(func)() = throw(MethodError($(func),()))
            function $(func)(::Int)
                q = Any[1]
                return q[1] # infer Any return type, workaround for julia issue #9479
            end
        end
        for t in tags
            t = quot(t)
            code = quote
                $code
                if haskey(methods_by_tag,$t)
                    push!(methods_by_tag[$t],$(quot(fname)))
                else
                    methods_by_tag[$t] = [$(quot(fname))]
                end
            end
        end
        push!(code.args, Expr(:export, func))
    end
    return code
end

abstract type AbstractMathProgModel end
export AbstractMathProgModel

# immutable type which we dispatch solvers on
abstract type AbstractMathProgSolver end
export AbstractMathProgSolver

# create dummy method to define function so that we can attach methods in other modules
# these are the common methods for AbstractMathProgModel
@define_interface begin
    getsolution
    getobjval
    optimize!
    status
    getobjbound, rewrap
    getobjgap, rewrap
    getrawsolver, rewrap
    getsolvetime, rewrap
    setsense!
    getsense
    numvar
    numconstr
    freemodel!
    setvartype!
    getvartype
    loadproblem!
    setwarmstart!
end

# solver parameters, may be implemented by AbstractMathProgModel or AbstractMathProgSolver
@define_interface begin
    setparameters!
end


include("LinearQuadratic.jl")
include("callbacks.jl")
include("Nonlinear.jl")
include("Conic.jl")

# Solver conversion routines
include("conic_to_lpqp.jl")
include("lpqp_to_conic.jl")
include("nonlinear_to_lpqp.jl")


end
