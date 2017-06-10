# Methods for the Nonlinear interface

abstract type AbstractNonlinearModel <: AbstractMathProgModel end
export AbstractNonlinearModel

abstract type AbstractNLPEvaluator end
export AbstractNLPEvaluator

# methods for AbstractNLPEvaluator
@define_interface begin
    NonlinearModel
    initialize
    features_available
    eval_f
    eval_g
    eval_grad_f
    jac_structure
    hesslag_structure
    eval_jac_g
    eval_jac_prod
    eval_jac_prod_t
    eval_hesslag_prod
    eval_hesslag
    isobjlinear
    isobjquadratic
    isconstrlinear
    obj_expr
    constr_expr
end

# fallback options
isobjlinear(::AbstractNLPEvaluator) = false
isobjquadratic(::AbstractNLPEvaluator) = false
isconstrlinear(::AbstractNLPEvaluator, i::Integer) = false

# getreducedcosts and getconstrduals already have stubs
