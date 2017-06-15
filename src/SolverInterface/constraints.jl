

"""
    ConstraintReference{T}

A lightweight object used to reference constraints in a model.
The parameter `T` is the type of set constraint referenced.
"""
struct ConstraintReference{T}
    value::UInt64
end

"""
    candelete(m::AbstractMathProgModel, ref::ConstraintReference)::Bool

Return a `Bool` indicating whether this constraint can be removed from the model `m`.
"""
candelete(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())

"""
    isvalid(m::AbstractMathProgModel, ref::ConstraintReference)::Bool

Return a `Bool` indicating whether this reference is valid
for an active constraint in the model `m`.
"""
isvalid(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())

"""
    delete!(m::AbstractMathProgModel, ref::ConstraintReference)

Delete the referenced constraint from the model.

    delete!(m::AbstractMathProgModel, refs::Vector{ConstraintReference})

Delete the referenced constraints in the vector `refs` from the model.
"""
Base.delete!(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())
Base.delete!(m::AbstractMathProgModel, refs::Vector{ConstraintReference}) = throw(MethodError())


"""
    addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}

Add the constraint
```math
Ax + b + q(x) \\in S
```
where ``A`` is a sparse vector specified in triplet form by
`a_constridx`, `a_varidx`, and `a_coef`; ``b`` is a vector;
``q(x)`` is a vector with component ``(q(x))_k`` defined to be ``\\frac{1}{2}x^TQ_kx``
where the symmetric matrix ``Q_k`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef` for which `Q_constridx` equals `k`; and the set ``S`` is defined by `S`.

Duplicate indices in either ``A`` or  a ``Q`` matrix are accepted and will be summed together. Off-diagonal entries of ``Q`` will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms. `a_varidx`, `Q_vari`, `Q_varj` should be collections of `VariableReference` objects.

    addconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}

A specialized version of `addconstraint!` for one-dimensional sets.
Add the constraint
```math
a^Tx + b + \\frac{1}{2}x^TQx \\in S
```
where ``a`` is a sparse vector specified in tuple form by
`a_varidx`, and `a_coef`; ``b`` is a scalar;
the symmetric matrix ``Q`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef`; and the set ``S`` is defined by `S`.

    addconstraint!(m::AbstractMathProgModel, varidx, S::AbstractSet)::ConstraintReference{typeof(S)}

A specialized version of `addconstraint!` for constraints on subsets of variables.
Add the constraint
```math
x_{varidx} \\in S
```
where `varidx` specifies the indices of the subvector of `x`.
"""
function addconstraint! end
