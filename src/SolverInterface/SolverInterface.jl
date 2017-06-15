using Base.Meta

export AbstractMathProgModel, AbstractMathProgSolver

# not sure if this should still be here
rewrap_methods = [:getobjbound,
                  :getobjgap,
                  :getrawsolver,
                  :getsolvetime,
                 ]

# immutable type which we dispatch solvers on
"""
    AbstractMathProgSolver

Abstract supertype for "solver" objects. A solver is a lightweight object used for selecting solvers and parameters. It does not store any instance data.
"""
abstract type AbstractMathProgSolver end

"""
    AbstractNLPModel

Abstract supertype which represents a solver's in-memory representation of a
non-linear optimization problem.
"""
abstract type AbstractNLPModel end

"""
    NLPModel(solver::AbstractMathProgSolver)

Create an instance of `AbstractNLPModel` using the given solver.
"""
function NLPModel end

"""
    AbstractModel

Abstract supertype which represents a solver's in-memory representation of an
optimization problem.
"""
abstract type AbstractModel end

"""
    Model(solver::AbstractMathProgSolver)

Create an instance of `AbstractModel` using the given solver.
"""
function Model end

"""
    AbstractMathProgModel

Union type of both `AbstractNLPModel` and `AbstractModel`.
"""
const AbstractMathProgModel = Union{AbstractNLPModel, AbstractModel}

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

"""
    writeproblem(m::AbstractModel, filename::String)

Writes the current problem data to the given file. Supported file types are solver-dependent.
"""
function writeproblem end

# these could change to attributes that can be set on a model *or solver*
# solver parameters, may be implemented by AbstractMathProgModel or AbstractMathProgSolver
function setparameters! end

include("variables.jl")
include("constraints.jl")
include("sets.jl")
include("attributes.jl")
include("objectives.jl")

#include("LinearQuadratic.jl")
#include("callbacks.jl")
#include("Nonlinear.jl")
#include("Conic.jl")

# Solver conversion routines
#include("conic_to_lpqp.jl")
#include("lpqp_to_conic.jl")
#include("nonlinear_to_lpqp.jl")
