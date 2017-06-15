```@meta
CurrentModule = MathProgBase
```

These are used to get and set properties of the model.

```@docs
AbstractAttribute
cangetattribute
getattribute
cansetattribute
setattribute!
```

## Scalar Attributes

```@docs
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
VariableCount
ConstraintCount
SupportsVariablesInSet
SupportsAffineInSet
SupportsQuadraticInSet
```

## Variable Attributes

These attributes are associated with variables. Calls to `getattribute` and `setattribute!` should include as an argument a single `VariableReference` or a vector of `VariableReference` objects.

```@docs
VariablePrimalStart
VariableLowerBoundDualStart
VariableUpperBoundDualStart
VariableLowerBound
VariableUpperBound
VariablePrimal
VariableLowerBoundDual
VariableUpperBoundDual
```

## Constraint Attributes

These attributes are associated with constraints. Calls to `getattribute` and `setattribute!` should include as an argument a single `ConstraintReference` or a vector of `ConstriaintReference{T}` objects.

```@docs
ConstraintPrimalStart
ConstraintDualStart
ConstraintPrimal
ConstraintDual
```


