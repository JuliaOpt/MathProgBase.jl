module SolverInterface

macro define_interface(args)
    @assert args.head == :block
    code = quote end
    for line in args.args
        isa(line, Symbol) || (line.head == :line && continue) || error("Unexpected code in block")
        func = esc(line)
        code = quote
            $code
            $(func)() = throw(MethodError($(func),()))
            function $(func)(::Int)
                q = Any[1]
                return q[1] # infer Any return type, workaround for julia issue #9479
            end
        end
        push!(code.args, Expr(:export, func))
    end
    return code
end

abstract AbstractMathProgModel
export AbstractMathProgModel

# immutable type which we dispatch solvers on 
abstract AbstractMathProgSolver
export AbstractMathProgSolver

# create dummy method to define function so that we can attach methods in other modules
@define_interface begin
    model
    loadproblem!
    writeproblem
    updatemodel!
    optimize!
end

include("setters.jl")
include("getters.jl")
include("callbacks.jl")
include("qcqp.jl")
include("sdp.jl")
include("nonlinear.jl")
include("conic.jl")
include("quad_to_conic.jl")

end
