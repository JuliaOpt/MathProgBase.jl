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
