```@meta
CurrentModule = MathProgBase
```

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
```
