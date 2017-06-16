```@meta
    CurrentModule = MathProgBase
```

# NonLinear Programming Interface (NLP)


## NLP Methods

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

## NLP Attributes

```@docs
ConstraintNLPDual
ConstraintNLPDualStart
```
