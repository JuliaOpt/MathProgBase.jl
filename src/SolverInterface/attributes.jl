"""
    AbstractAttribute

Abstract supertype for attribute objects that can be used
to set or get attributes (properties) of the model.
"""
abstract type AbstractAttribute end

"""
    getattribute(m::AbstractMathProgModel, attr::AbstractAttribute, extra_args...)

Return an attribute of the model `m` specified by attribute type `attr`.

# Examples

    getattribute(m, ObjectiveValue())
    getattribute(m, VariableResult(), ref)
    getattribute(m, VariableResult(5), [ref1,ref2])
    getattribute(m, OtherAttribute("something specific to cplex"))
"""
function getattribute end

function getattribute(m, attr::AbstractAttribute, args...)
    throw(ArgumentError("Model of type $(typeof(m)) does not support accessing the attribute $attr"))
end

"""
    getattribute!(output, m::AbstractMathProgModel, attr::AbstractAttribute, extra_args...)

An in-place version of `getattribute`. Return an attribute of the model `m` specified by attribute type `attr` into output vector `attr`.
"""
function getattribute! end
function getattribute!(output, m, attr::AbstractAttribute, args...)
    throw(ArgumentError("Model of type $(typeof(m)) does not support accessing the attribute $attr"))
end

"""
    cangetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool

Return a `Bool` indicating whether the model `m` currently has a value for
the attributed specified by attribute type `attr`.
"""
function cangetattribute end
cangetattribute(m::AbstractMathProgModel, attr::AbstractAttribute) = false

"""
    cansetattribute(m::AbstractMathProgModel, attr::AbstractAttribute)::Bool

Return a `Bool` indicating whether the model `m` will accept a
`setattribute!` call for the attributed specified by attribute type `attr`.
"""
function cansetattribute end
cansetattribute(m::AbstractMathProgModel, attr::AbstractAttribute) = false

"""
    setattribute!(m::AbstractMathProgModel, attr::AbstractAttribute, ...)

Set an attribute of the model `m` specified by attribute type `attr`.
"""
function setattribute! end
function setattribute!(m, attr::AbstractAttribute, args...)
    throw(ArgumentError("Model of type $(typeof(m)) does not support setting the attribute $attr"))
end

# Scalar attributes
"""
    ObjectiveValue(N)
    ObjectiveValue()

The objective value of the primal result `N`. If `N` is omitted, it defaults to 1.
"""
struct ObjectiveValue <: AbstractAttribute
    N::Int
end
ObjectiveValue() = ObjectiveValue(1)

"""
    ObjectiveBound()

The best known bound on the optimal objective value.
"""
struct ObjectiveBound <: AbstractAttribute end

"""
    RelativeGap()

The final relative optimality gap as optimization terminated. That is, ``\\frac{|b-f|}{|f|}``, where ``b`` is the best bound and ``f`` is the best feasible objective value.
"""
struct RelativeGap <: AbstractAttribute  end

"""
    SolveTime()

The total elapsed solution time (in seconds) as reported by the solver.
"""
struct SolveTime <: AbstractAttribute end

"""
    Sense()

The optimization sense of the model, an `OptimizationSense` with value `MinSense` or `MaxSense`.
"""
struct Sense <: AbstractAttribute end

@enum OptimizationSense MinSense MaxSense

"""
    SimplexIterations()

The cumulative number of simplex iterations during the optimization process. In particular, for a MIP the total simplex iterations for all nodes.
"""
struct SimplexIterations <: AbstractAttribute end

"""
    BarrierIterations()

The cumulative number of barrier iterations during the optimization process.
"""
struct BarrierIterations <: AbstractAttribute end

"""
    NodeCount()

The total number of branch-and-bound nodes explored.
"""
struct NodeCount <: AbstractAttribute end

"""
    RawSolver()

An object that may be used to access a solver-specific API for this model.
"""
struct RawSolver <: AbstractAttribute end

"""
    ResultCount()

The number of results available.
"""
struct ResultCount <: AbstractAttribute end

"""
    VariableCount()

The number of variables in the model.
"""
struct VariableCount <: AbstractAttribute end

"""
    ConstraintCount{T}()

The number of constraints of type T in the model.
"""
struct ConstraintCount{T} <: AbstractAttribute end

"""
    SupportsVariablesInSet{T}()

A `Bool` indicating whether the solver or model supports
a constraint of type ``x_{varidx} \\in S`` where ``S``
is a set of type `T` and ``varidx`` indicates any subset
of the variables.
"""
struct SupportsVariablesInSet{T} <: AbstractAttribute end

"""
    SupportsAffineInSet{T}()

A `Bool` indicating whether the solver or model supports
a constraint of of the form "affine expression" in ``S`` where ``S``
is a set of type `T`.
"""
struct SupportsAffineInSet{T} <: AbstractAttribute end

"""
    SupportsQuadraticInSet{T}()

A `Bool` indicating whether the solver or model supports
a constraint of of the form "quadratic expression" in ``S`` where ``S``
is a set of type `T`.
"""
struct SupportsQuadraticInSet{T} <: AbstractAttribute end



# Variable attributes

"""
    VariablePrimalStart()

An initial assignment of the variables that the solver may use
to warm-start the solve.
"""
struct VariablePrimalStart <: AbstractAttribute end

"""
    VariableLowerBoundDualStart()

An initial assignment of the variable lower-bound duals that the solver may use
to warm-start the solve.
"""
struct VariableLowerBoundDualStart <: AbstractAttribute end

"""
    VariableUpperBoundDualStart()

An initial assignment of the variable upper-bound duals that the solver may use
to warm-start the solve.
"""
struct VariableUpperBoundDualStart <: AbstractAttribute end

"""
    VariableLowerBound()

Lower-bound constraints on variables. `-Inf` is valid as no bound.
"""
struct VariableLowerBound <: AbstractAttribute end

"""
    VariableUpperBound()

Upper-bound constraints for the variables. `Inf` is valid as no bound.
"""
struct VariableUpperBound <: AbstractAttribute end

"""
    VariablePrimal(N)
    VariablePrimal()

The assignment to the primal variables in result `N`. If `N` is omitted, it is 1 by default.
"""
struct VariablePrimal <: AbstractAttribute
    N::Int
end
VariablePrimal() = VariablePrimal(1)

"""
    VariableLowerBoundDual(N)

The assignment to the duals on the variable lower bounds in result `N`. If `N` is omitted, it is interpreted as 1.
"""
struct VariableLowerBoundDual <: AbstractAttribute
    N::Int
end
VariableLowerBoundDual() = VariableLowerBoundDual(1)

"""
    VariableUpperBoundDual(N)

The assignment to the duals on the variable upper bounds in result `N`. If `N` is omitted, it is interpreted as 1.
"""
struct VariableUpperBoundDual <: AbstractAttribute
    N::Int
end
VariableUpperBoundDual() = VariableUpperBoundDual(1)

"""
    VariableBasisStatus()

Returns the `VariableBasisStatusCode` of a given variable, with respect to an available optimal solution basis.
    * `Basic`: variable is in the basis.
    * `NonbasicAtLower`: variable is not in the basis and is at its lower bound.
    * `NonbasicAtUpper`: variable is not in the basis and is at its upper bound.
    * `SuperBasic`: variable is not in the basis but is also not at one of its bounds.
"""
struct VariableBasisStatus <: AbstractAttribute end
@enum VariableBasisStatusCode Basic NonbasicAtLower NonbasicAtUpper SuperBasic

# Constraint attributes

"""
    ConstraintPrimalStart()

An initial assignment of the constraint primal values that the solver may use
to warm-start the solve.
"""
struct ConstraintPrimalStart <: AbstractAttribute end

"""
    ConstraintDualStart()

An initial assignment of the constriant duals that the solver may use
to warm-start the solve.
"""
struct ConstraintDualStart <: AbstractAttribute end

"""
    ConstraintPrimal(N)
    ConstraintPrimal()

The assignment to the constraint primal values in result `N`. If `N` is omitted, it is 1 by default.
"""
struct ConstraintPrimal <: AbstractAttribute
    N::Int
end
ConstraintPrimal() = ConstraintPrimal(1)

"""
    ConstraintDual(N)
    ConstraintDual()

The assignment to the constraint dual values in result `N`. If `N` is omitted, it is 1 by default.
"""
struct ConstraintDual <: AbstractAttribute
    N::Int
end
ConstraintDual() = ConstraintDual(1)

"""
    ConstraintBasisStatus()

Returns the `ConstraintBasisStatusCode` of a given constraint, with respect to an available optimal solution basis.
    * `Basic`: constraint is in the basis, i.e. the corresponding slack variable is in the basis.
    * `Nonbasic`: constraint is not in the basis.
"""
struct ConstraintBasisStatus <: AbstractAttribute end
@enum ConstraintBasisStatusCode Basic Nonbasic




# Termination status
"""
    TerminationStatus()

A `TerminationStatusCode` explaining why the solver stopped.
"""
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

# Result status

"""
    ResultStatusCode

An Enum of possible values for the `PrimalStatus` and `DualStatus` attributes. The values indicate how to interpret the result vector.

  * `FeasiblePoint`
  * `NearlyFeasiblePoint`
  * `InfeasiblePoint`
  * `InfeasibilityCertificate`
  * `NearlyInfeasibilityCertificate`
  * `ReductionCertificate`
  * `NearlyReductionCertificate`
  * `Unknown`
  * `Other`
"""
@enum ResultStatusCode FeasiblePoint NearlyFeasiblePoint InfeasiblePoint InfeasibilityCertificate NearlyInfeasibilityCertificate ReductionCertificate NearlyReductionCertificate Unknown Other

"""
    PrimalStatus(N)
    PrimalStatus()

The `ResultStatusCode` of the primal result `N`. If `N` is omitted, it defaults to 1.
"""
struct PrimalStatus <: AbstractAttribute
    N::Int
end
PrimalStatus() = PrimalStatus(1)

"""
    DualStatus(N)
    DualStatus()

The `ResultStatusCode` of the dual result `N`. If `N` is omitted, it defaults to 1.
"""
struct DualStatus <: AbstractAttribute
    N::Int
end
DualStatus() = DualStatus(1)
