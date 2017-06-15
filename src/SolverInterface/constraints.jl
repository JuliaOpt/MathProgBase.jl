




"""
    addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}

Add the quadratic-in-set constraint
```math
Ax + b + q(x) \\in S
```
where ``A`` is a sparse matrix specified in triplet form by
`a_constridx`, `a_varidx`, and `a_coef`; ``b`` is a vector;
``q(x)`` is a vector with component ``(q(x))_k`` defined to be ``\\frac{1}{2}x^TQ_kx``
where the symmetric matrix ``Q_k`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef` for which `Q_constridx` equals `k`; and the set ``S`` is defined by `S`.

Duplicate indices in either ``A`` or  a ``Q`` matrix are accepted and will be summed together. Off-diagonal entries of ``Q`` will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms. `a_varidx`, `Q_vari`, `Q_varj` should be collections of `VariableReference` objects.



    addconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)::QuadraticConstraintReference{typeof(S)}

A specialized version of `addconstraint!` for one-dimensional sets.
Add the constraint
```math
a^Tx + b + \\frac{1}{2}x^TQx \\in S
```
where ``a`` is a sparse vector specified in tuple form by
`a_varidx`, and `a_coef`; ``b`` is a scalar;
the symmetric matrix ``Q`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef`; and the set ``S`` is defined by `S`.

    addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}

Add the affine-in-set constraint
```math
Ax + b \\in S
```
where ``A`` is a sparse matrix specified in triplet form by
`a_constridx`, `a_varidx`, and `a_coef`; ``b`` is a vector; and the set ``S`` is defined by `S`.

Duplicate indices either ``A`` are accepted and will be summed together.

    addconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, S::AbstractSet)::AffineConstraintReference{typeof(S)}

A specialized version of `addconstraint!` for one-dimensional sets.
Add the constraint
```math
a^Tx + b \\in S
```
where ``a`` is a sparse vector specified in tuple form by
`a_varidx`, and `a_coef`; ``b`` is a scalar; and the set ``S`` is defined by `S`.


    addconstraint!(m::AbstractMathProgModel, varidx, S::AbstractSet)::VariablewiseConstraintReference{typeof(S)}

A specialized version of `addconstraint!` for variablewise constraints.
Add the constraint
```math
x_{varidx} \\in S
```
where `varidx` specifies the indices of the subvector of `x`.
"""
function addconstraint! end

"""
    modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, args...)

Modify elements of the `i`'th row of the constraint `c` depending on the
arguments `args`. The `i`'th row will have the form
```math
    a_i^T + b_i + \\frac{1}{2}x^TQ_ix \\in S
```
There are three cases.

# Modify Constant term

    modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, b)

Set the constant term of the `i`'th row in the constraint `c` to `b`.

### Examples

        modifyconstraint!(m, c, 1, 1.0)

# Modify Linear term

modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, a_varidx, a_coef)

Set elements given by `a_varidx` in the linear term of the `i`'th element in the
constraint `c` to `a_coef`. Either `a_varidx` and `a_coef` are both singletons,
or they should be collections with equal length.

The behaviour of duplicate entries in `a_varidx` is undefined.

### Examples

        modifyconstraint!(m, c, v, 1.0)
        modifyconstraint!(m, c, [v1, v2], [1.0, 2.0])

# Modify Quadratic term

    modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference, i::Int, Q_vari, Q_varj, Q_coef)

Set the elements in the quadratic term of the `i`'th element of the constraint `c`
specified by the triplets `Q_vari`, `Q_varj`, and `Q_coef`. Off-diagonal entries
will be mirrored. `Q_vari`, `Q_varj` should be collections of `VariableReference`
objects.

The behaviour of duplicate entries is undefined. If entries for both ``(i,j)``
and ``(j,i)`` are provided, these are considered duplicate terms.

### Examples

        modifyconstraint!(m, c, v1, v2, 1.0)
        modifyconstraint!(m, c, [v1, v2], [v1, v1], [1.0, 2.0])

# Modify Set

    modifyconstraint!(m::AbstractMathProgModel, c::ConstraintReference{S}, set::S)

Change the set of constraint `c` to the new set `set` which should be of the same
type as the original set.

### Examples

If `c` is a `ConstraintReference{Interval}`

        modifyconstraint!(m, c, Interval(0, 5))
        modifyconstraint!(m, c, NonPositive) # errors
"""
function modifyconstraint! end
