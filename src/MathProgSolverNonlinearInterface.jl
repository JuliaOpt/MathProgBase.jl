

export AbstractNLPEvaluator


abstract AbstractNLPEvaluator

# methods for AbstractNLPEvaluator
for func in [:initialize,
             :features_available,
             :eval_f,
             :eval_g,
             :eval_grad_f,
             :jac_structure,
             :hesslag_structure,
             :eval_jac_g,
             :eval_jac_prod,
             :eval_jac_prod_t,
             :eval_hesslag_prod,
             :eval_hesslag,
             :isobjlinear,
             :isobjquadratic,
             :isconstrlinear,
             ] # :obj_expr and :constr_expr not yet defined
    @eval $(func)() = throw(MethodError($(func),()))
    eval(Expr(:export,func))
end

# additional methods for AbstractMathProgModel
loadnonlinearproblem!() = throw(MethodError(:loadnonlinearproblem!,()))
export loadnonlinearproblem!
