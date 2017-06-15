"""
    setobjective!(m::AbstractMathProgModel, N::Int, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef)

Set the `N`'th objective in the model `m` to be
```math
a^Tx + b + \\frac{1}{2}x^TQx
```
where ``a`` is a sparse vector specified in tuple form by `a_varidx`, and
`a_coef`; ``b`` is a scalar; and the symmetric matrix ``Q`` is defined by the
triplets in `Q_vari`, `Q_varj`, `Q_coef`.

Duplicate indices in either ``A`` or  a ``Q`` matrix are accepted and will be
summed together. Off-diagonal entries of ``Q`` will be mirrored, so either the
upper triangular or lower triangular entries of ``Q`` should be provided. If
entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered
duplicate terms. `a_varidx`, `Q_vari`, `Q_varj` should be collections of
`VariableReference` objects.
"""
function setobjective! end

"""
    modifyobjective!(m::AbstractMathProgModel, i::Int, args...)

Modify elements of the `i`'th objective depending on the
arguments `args`. The `i`'th objective will have the form:
```math
    a_i^T + b_i + \\frac{1}{2}x^TQ_ix
```
There are three cases.

# Modify Constant term

    modifyobjective!(m::AbstractMathProgModel, i::Int, b)

Set the constant term of the `i`'th row objective to `b`.

### Examples

        modifyobjective!(m, 1, 1.0)

# Modify Linear term

    modifyobjective!(m::AbstractMathProgModel, i::Int, a_varidx, a_coef)

Set elements given by `a_varidx` in the linear term of the `i`'th objective to
`a_coef`. Either `a_varidx` and `a_coef` are both singletons, or they should be
collections with equal length.

The behaviour of duplicate entries in `a_varidx` is undefined.

### Examples

        modifyobjective!(m, 1, v, 1.0)
        modifyobjective!(m, 1, [v1, v2], [1.0, 2.0])

# Modify Quadratic term

    modifyobjective!(m::AbstractMathProgModel, i::Int, Q_vari, Q_varj, Q_coef)

Set the elements in the quadratic term of the `i`'th objective specified by the
triplets `Q_vari`, `Q_varj`, and `Q_coef`. Off-diagonal entries will be mirrored.
`Q_vari`, `Q_varj` should be collections of `VariableReference` objects.

The behaviour of duplicate entries is undefined. If entries for both ``(i,j)``
and ``(j,i)`` are provided, these are considered duplicate terms.

### Examples

        modifyobjective!(m, 1, v1, v2, 1.0)
        modifyobjective!(m, 1, [v1, v2], [v1, v1], [1.0, 2.0])
"""
function modifyobjective! end

"""
    getobjective(m, i:Int)

Returns the `i`'th objective as the tuple `(b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef)`.

The elements in the tuple are the same as those defined in `addobjective!`.
"""
function getobjective end
