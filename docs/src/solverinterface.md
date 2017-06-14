# Solver Interface

```@meta
CurrentModule = MathProgBase
```

```@docs
AbstractMathProgModel
AbstractMathProgSolver
```

## Basic Methods

```@docs
optimize!
freemodel!
```

## Sets

How to add constraints.
```@docs
addconstraint!
```

List of sets.
```@docs
NonNegative
NonPositive
Zero
Interval
```

## Attributes

These are used to get and set properties of the model.

```@docs
AbstractAttribute
cangetattribute
getattribute
cansetattribute
setattribute!
```

### Scalar Attributes

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
```

### Variable Attributes

```@docs
VariableStart
VariableDualStart
VariableResult
```


## Termination Status

The `TerminationStatus` attribute is meant to explain the reason why the solver stopped executing. The value of the attribute is of type `TerminationStatusCode`.

```@docs
TerminationStatusCode
```
