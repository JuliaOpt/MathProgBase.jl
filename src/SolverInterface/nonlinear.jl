abstract AbstractNLPEvaluator
export AbstractNLPEvaluator

# methods for AbstractNLPEvaluator
@define_interface begin
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
end

# additional methods for AbstractMathProgModel
loadnonlinearproblem!() = throw(MethodError(:loadnonlinearproblem!,()))
export loadnonlinearproblem!
