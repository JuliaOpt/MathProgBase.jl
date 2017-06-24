```@meta
CurrentModule = MathProgBase
```

# Solver Interface API

Some introduction to MPB API. List basic standalone methods.

```@docs
AbstractModel
AbstractNLPModel
AbstractMathProgSolver
```


```@docs
Model
NLPModel
optimize!
freemodel!
```

## Variables


```@docs
VariableReference
candelete(::AbstractMathProgModel,::VariableReference)
isvalid(::AbstractMathProgModel,::VariableReference)
delete!(::AbstractMathProgModel,::VariableReference)
addvariables!
addvariable!
```

## Objectives

How to add and set objectives.
```@docs
setobjective!
modifyobjective!
getobjectiveconstant
getobjectiveaffine
```

## Constraints

How to add and modify constraints.
```@docs
VariablewiseConstraintReference
AffineConstraintReference
QuadraticConstraintReference
candelete(::AbstractMathProgModel,::ConstraintReference)
isvalid(::AbstractMathProgModel,::ConstraintReference)
delete!(::AbstractMathProgModel,::ConstraintReference)
addconstraint!
modifyconstraint!
getconstraintconstant
getconstraintaffine
getconstraintquadratic
```

## Sets

List of sets.
```@docs
NonNegative
NonPositive
Zero
Interval
SecondOrderCone
PositiveSemidefiniteConeTriangle
PositiveSemidefiniteConeScaled
Integers
Binaries
SOS1
SOS2
```

## Attributes

These are used to get and set properties of the model or solver.

```@docs
AbstractSolverOrModelAttribute
AbstractVariableAttribute
AbstractConstraintAttribute
cangetattribute
getattribute
getattribute!
cansetattribute
setattribute!
```

### Solver or Model Attributes

These attribute apply to solver or model objects, as noted.

```@docs
ReturnsDuals
SupportsAddConstraintAfterSolve
SupportsDeleteConstraint
SupportsAddVariableAfterSolver
SupportsQuadraticObjective
SupportsConicThroughQuadratic
ObjectiveValue
ObjectiveBound
RelativeGap
SolveTime
Sense
SimplexIterations
BarrierIterations
NodeCount
RawSolver
ResultCount
NumberOfVariables
NumberOfVariablewiseConstraints
NumberOfAffineConstraints
NumberOfQuadraticConstraints
SupportsVariablewiseConstraint
SupportsAffineConstraint
SupportsQuadraticConstraint
TerminationStatus
PrimalStatus
DualStatus
```

### Variable Attributes

These attributes are associated with variables. Calls to `getattribute` and `setattribute!` should include as an argument a single `VariableReference` or a vector of `VariableReference` objects. By convention, the name of a variable attributes begins with `Variable`.

```@docs
VariablePrimalStart
VariableLowerBoundDualStart
VariableUpperBoundDualStart
VariableLowerBound
VariableUpperBound
VariablePrimal
VariableLowerBoundDual
VariableUpperBoundDual
VariableBasisStatus
```

### Constraint Attributes

These attributes are associated with constraints. Calls to `getattribute` and `setattribute!` should include as an argument a single `ConstraintReference` or a vector of `ConstriaintReference{T}` objects. By convention, the name of a variable attributes begins with `Variable`.


```@docs
ConstraintPrimalStart
ConstraintDualStart
ConstraintPrimal
ConstraintDual
ConstraintBasisStatus
```

## Status Codes

### Termination Status

The `TerminationStatus` attribute is meant to explain the reason why the solver stopped executing. The value of the attribute is of type `TerminationStatusCode`.

```@docs
TerminationStatusCode
```

### Result Status

The `PrimalStatus` and `DualStatus` attributes are meant to explain how to interpret the result returned by the solver. The value of the attributes are of type `ResultStatusCode`.

```@docs
ResultStatusCode
```

### Basis Status

```@docs
BasisStatusCode
```

## Duals

We take the convention that duals on variable lower bounds should be nonnegative, duals on variable upper bounds should be nonpositive, and duals on closed convex cones should belong to the dual cone.

## Nonlinear Programming (NLP)


### NLP Methods

```@docs
loadnlp!
initialize
features_available
eval_f
eval_grad_f
jac_structure
hesslag_structure
eval_jac_g
eval_jac_prod
eval_jac_prod_t
eval_hesslag_prod
eval_hesslag
isobjlinear(::AbstractNLPEvaluator)
isobjquadratic(::AbstractNLPEvaluator)
isconstrlinear(::AbstractNLPEvaluator, i::Integer)
obj_expr
constr_expr
```

### NLP Attributes

```@docs
ConstraintNLPDual
ConstraintNLPDualStart
```
