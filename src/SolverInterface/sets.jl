

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

