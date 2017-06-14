"""
    AbstractAttribute

Abstract supertype for attribute objects that can be used
to set or get attributes (properties) of the model.
"""
abstract type AbstractAttribute end

# TODO:
# in-place getters

"""
    getattribute(m::AbstractMathProgModel, attr::AbstractAttribute)

Return an attribute of the model `m` specified by attribute type `attr`.

# Examples

    getattribute(m, ObjectiveValue())
    getattribute(m, VariableResult())
    getattribute(m, VariableResult(5))
    getattribute(m, OtherAttribute("something specific to cplex"))
"""
function getattribute end

"""
    cangetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool

Return a `Bool` indicating whether the model `m` currently has a value for
the attributed specified by attribute type `attr`.
"""
function cangetattribute end

"""
    cansetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool

Return a `Bool` indicating whether the model `m` will accept a
`setattribute!` call for the attributed specified by attribute type `attr`.
"""
function cansetattribute end

"""
    setattribute!(m::AbstractMathProgModel, attr::AbstractAttribute, ...)

Set an attribute of the model `m` specified by attribute type `attr`.
"""
function setattribute! end

# Scalar attributes
"""
    ObjectiveValue

The objective value of the best primal result vector found by the solver.
"""
struct ObjectiveValue <: AbstractAttribute end

"""
    ObjectiveBound

The best known bound on the optimal objective value.
"""
struct ObjectiveBound <: AbstractAttribute end


struct RelativeGapAttr <: AbstractAttribute  end
"""
    RelativeGap

The final relative optimality gap as optimization terminated. That is, ``\\frac{|b-f|}{|f|}``, where ``b`` is the best bound and ``f`` is the best feasible objective value.
"""
const RelativeGap = RelativeGapAttr()

"""
    SolveTime

The total elapsed solution time (in seconds) as reported by the solver.
"""
struct SolveTime <: AbstractAttribute end

"""
    Sense

The optimization sense of the model. Valid values are `:Min` and `:Max`.
"""
struct Sense <: AbstractAttribute end

"""
    SimplexIterations

The cumulative number of simplex iterations during the optimization process. In particular, for a MIP the total simplex iterations for all nodes.
"""
struct SimplexIterations <: AbstractAttribute end

"""
    BarrierIterations

The cumulative number of barrier iterations during the optimization process.
"""
struct BarrierIterations <: AbstractAttribute end

"""
    NodeCount

The total number of branch-and-bound nodes explored.
"""
struct NodeCount <: AbstractAttribute end

"""
    RawSolver

An object that may be used to access a solver-specific API for this model.
"""
struct RawSolver <: AbstractAttribute end

"""
    ResultCount

The number of results available.
"""
struct ResultCount <: AbstractAttribute end


# Variable attributes

"""
    VariableStart

An initial assignment of the variables that the solver may use
to warm-start the solve.
"""
struct VariableStart <: AbstractAttribute end

"""
    VariableDualStart

An initial assignment of the variable duals that the solver may use
to warm-start the solve.
"""
struct VariableDualStart <: AbstractAttribute end

"""
    VariableResult

The assignment to the primal variables in result `N`. If `N` is omitted, it is interpreted as 1.
"""
struct VariableResult <: AbstractAttribute
    N::Int
end

# VarType?


function getsolution end
function getobjval end


function setvartype! end
function getvartype end
function loadproblem! end
function setwarmstart! end


# Termination status

struct TerminationStatus <: AbstractAttribute end

"""
    TerminationStatusCode

An Enum of possible values for the `TerminationStatus` attribute. This attribute is meant to explain the reason why the solver stopped executing.

# OK

These are generally OK statuses.

  * `Success`: the algorithm ran successfully and has a result. This includes cases where the algorithm converges to an infeasible point (NLP) or converges to a solution of a homogeneous self-dual problem and has a certificate of primal/dual infeasibility.

  * `AlmostSuccess`: the algorithm *almost* ran successfully (e.g., to relaxed convergence tolerances) and has a result.

  * `InfeasibleNoResult`: the algorithm stopped because it decided that the problem is infeasible but does not have a result to return.

  * `UnboundedNoResult`: the algorithm stopped because it decided that the problem is unbounded but does not have a result to return.

  * `InfeasibleOrUnbounded`: the algorithm stopped because it decided that the problem is infeasible or unbounded; no result is available. This occasionally happens during MIP presolve.

# Limits

The solver stopped because of some user-defined limit.
To be documented: `IterationLimit`, `TimeLimit`, `NodeLimit`, `SolutionLimit`, `MemoryLimit`, `ObjectiveLimit`, `NormLimit`, `OtherLimit`.

# Problematic

This group of statuses means that something unexpected or problematic happened.

  * `SlowProgress`: the algorithm stopped because it was unable to continue making progress towards the solution. `AlmostSuccess` should be used if there is additional information that relaxed convergence tolerances are satisfied.

To be documented: `NumericalError`, `InvalidModel`, `InvalidOption`, `Interrupted`, `OtherError`.

"""
@enum TerminationStatusCode Success AlmostSuccess InfeasibleNoResult UnboundedNoResult InfeasibleOrUnbounded IterationLimit TimeLimit NodeLimit SolutionLimit MemoryLimit ObjectiveLimit NormLimit OtherLimit SlowProgress NumericalError InvalidModel InvalidOption Interrupted OtherError


