

"""
    AbstractSet
    
Abstract supertype for set objects used to encode constraints.
"""
abstract type AbstractSet end


"""
    NonNegative(n)

The nonnegative orthant ``\\{ x \\in \\mathbb{R}^n : x \\ge 0 \\}`` where the dimension ``n`` is specified by the field `n`.
"""
struct NonNegative <: AbstractSet
    dim::Int
end

"""
    NonPositive(n)

The nonpositive orthant ``\\{ x \\in \\mathbb{R}^n : x \\le 0 \\}`` where the dimension ``n`` is specified by the field `n`.
"""
struct NonPositive <: AbstractSet
    dim::Int
end

"""
    Zero(n)

The set ``\\{0\\}^n`` where the dimension ``n`` is specified by the field `n`.
"""
struct Zero <: AbstractSet
    dim::Int
end


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
    Integers(n)

The set of integers ``\\mathbb{Z}^n``.
"""
struct Integers <: AbstractSet
    dim::Int
end

"""
    Binaries(n)

The set of binary vectors ``\\{0,1\\}^n``.
"""
struct Binaries <: AbstractSet
    dim::Int
end


dimension(s::Union{NonNegative,NonPositive,Zero,Integers,Binaries}) = s.dim

"""
    addconstraint!(m::AbstractMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)

Add the constraint
```math
Ax + b + q(x) \\in S
```
where ``A`` is a sparse vector specified in triplet form by
`a_constridx`, `a_varidx`, and `a_coef`; ``b`` is a vector;
``q(x)`` is a vector with component ``(q(x))_k`` defined to be ``\\frac{1}{2}x^TQ_kx``
where the symmetric matrix ``Q_k`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef` for which `Q_constridx` equals `k`; and the set ``S`` is defined by `S`.

Duplicate indices in either ``A`` or  a ``Q`` matrix are accepted and will be summed together. Off-diagonal entries of ``Q`` will be mirrored, so either the upper triangular or lower triangular entries of ``Q`` should be provided. If entries for both ``(i,j)`` and ``(j,i)`` are provided, these are considered duplicate terms.

    addconstraint!(m::AbstractMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)

A specialized version of `addconstraint!` for one-dimensional sets.
Add the constraint
```math
a^Tx + b + \\frac{1}{2}x^TQx \\in S
```
where ``a`` is a sparse vector specified in tuple form by
`a_varidx`, and `a_coef`; ``b`` is a scalar;
the symmetric matrix ``Q`` is defined by the triplets in `Q_vari`, `Q_varj`,
`Q_coef`; and the set ``S`` is defined by `S`.

    addconstraint!(m::AbstractMathProgModel, varidx, S::AbstractSet)

A specialized version of `addconstraint!` for constraints on subsets of variables.
Add the constraint
```math
x_{varidx} \\in S
```
where `varidx` specifies the indices of the subvector of `x`.
"""
function addconstraint! end
