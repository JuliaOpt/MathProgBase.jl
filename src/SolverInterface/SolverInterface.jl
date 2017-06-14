using Base.Meta

export AbstractMathProgModel, AbstractMathProgSolver

# not sure if this should still be here
rewrap_methods = [:getobjbound,
                  :getobjgap,
                  :getrawsolver,
                  :getsolvetime,
                 ]

"""
    AbstractMathProgModel

Abstract supertype which represents a solver's
in-memory representation of an optimization problem.
"""
abstract type AbstractMathProgModel end

# immutable type which we dispatch solvers on
"""
    AbstractMathProgSolver

Abstract supertype for "solver" objects. A solver is a lightweight object used for selecting solvers and parameters. It does not store any instance data.
"""
abstract type AbstractMathProgSolver end


# basic methods methods for AbstractMathProgModel

"""
    optimize!(m::AbstractMathProgModel)

Start the solution procedure.
"""
function optimize! end

"""
    freemodel!(m::AbstractMathProgModel)

Release any resources and memory used by the model. Note that the
Julia garbage collector takes care of this automatically, but
automatic collection cannot always be forced. This method is useful for more
precise control of resources, especially in the case of commercial solvers
with licensing restrictions on the number of concurrent runs.
Users must discard the model object after this method is invoked.
"""
function freemodel! end

# TODO
function loadproblem! end



# these could change to attributes that can be set on a model *or solver*
# solver parameters, may be implemented by AbstractMathProgModel or AbstractMathProgSolver
function setparameters! end

include("variables.jl")
include("constraints.jl")
include("sets.jl")
include("attributes.jl")

#include("LinearQuadratic.jl")
#include("callbacks.jl")
#include("Nonlinear.jl")
#include("Conic.jl")

# Solver conversion routines
#include("conic_to_lpqp.jl")
#include("lpqp_to_conic.jl")
#include("nonlinear_to_lpqp.jl")
