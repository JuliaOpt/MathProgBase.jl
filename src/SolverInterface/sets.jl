

"""
    AbstractSet
    
Abstract supertype for set objects used to encode constraints.
"""
abstract type AbstractSet end


"""
    NonNegative(dim)

The nonnegative orthant ``\\{ x \\in \\mathbb{R}^n : x \\ge 0 \\}`` where the dimension ``n`` is specified by the field `dim`.
"""
struct NonNegative <: AbstractSet
    dim::Int
end

"""
    NonPositive(dim)

The nonpositive orthant ``\\{ x \\in \\mathbb{R}^n : x \\le 0 \\}`` where the dimension ``n`` is specified by the field `dim`.
"""
struct NonPositive <: AbstractSet
    dim::Int
end

"""
    Zero(dim)

The set ``\\{0\\}^n`` where the dimension ``n`` is specified by the field `dim`.
"""
struct Zero <: AbstractSet
    dim::Int
end

dimension(s::Union{NonNegative,NonPositive,Zero}) = s.dim

"""
    Interval(lower,upper)

The set ``[l,u] \\subseteq \\mathbb{R}^n`` where ``l`` and ``u`` are specified by `lower` and `upper`, respectively. We allow `lower` and `upper` to be `-Inf` or `Inf`, in which case the set is interpreted as a one-sided interval.
"""
struct Interval{T} <: AbstractSet
    lower::T
    upper::T
end

dimension(s::Interval) = 1

"""
    addconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)

Add the constraint
```math
a^Tx + b + \\frac{1}{2}x^TQx \\in S
```
where ``a`` is a sparse vector specified in tuple form by
`a_varidx`, and `a_coef`; ``b`` is a scalar;
the symmetric matrix ``Q_i`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef`; and the set ``S`` is defined by `S`.

Duplicate indices in either `a_varidx` or `Q_vari` and `Q_varj` are accepted and will be summed together. Off-diagonal entries of ``Q`` will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms.
"""
function addconstraint! end
