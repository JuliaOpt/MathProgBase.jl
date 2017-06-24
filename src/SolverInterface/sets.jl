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

dimension(s::Interval) = length(s.lower)

"""
    SecondOrderCone(n)

The second-order cone or the Lorenz cone of dimension `n`
defined as
```math
\\{ (t,x) \\in \\mathbb{R}^n : t \\ge ||x||_2 \\}.
```
"""
struct SecondOrderCone <: AbstractSet
    dim::Int
end

#ExponentialCone

"""
    PositiveSemidefiniteConeTriangle(n)

The cone of symmetric ``n \\times n`` matrices that are positive semidefinite.
The dimension of the cone is ``n(n+1)/2`` since the matrices are symmetric.
The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column).
The scalar product is the sum of the pairwise product of the diagonal entries plus twice the sum of the pairwise product of the upper diagonal entries.

### Examples

The matrix
```math
\\begin{bmatrix}
  1 & 2 & 3\\\\
  2 & 4 & 5\\\\
  3 & 5 & 6
\\end{bmatrix}
```
corresponds to ``(1, 2, 3, 4, 5, 6)`` for `PositiveSemidefiniteConeTriangle`
"""
struct PositiveSemidefiniteConeTriangle <: AbstractSet
    dim::Int
end

"""
    PositiveSemidefiniteConeScaled(n)

The cone of symmetric ``n \\times n`` matrices that are positive semidefinite.
The dimension of the cone is ``n(n+1)/2`` since the matrices are symmetric.
The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column).
The off-diagonal entries of the matrices of both the cone and its dual are scaled by ``\\sqrt{2}`` and the scalar product is simply the sum of the pairwise product of the entries.

### Examples

The matrix
```math
\\begin{bmatrix}
  1 & 2 & 3\\\\
  2 & 4 & 5\\\\
  3 & 5 & 6
\\end{bmatrix}
```
and to ``(1, 2\\sqrt{2}, 3\\sqrt{2}, 4, 5\\sqrt{2}, 6)`` for `PositiveSemidefiniteConeScaled`.
"""
struct PositiveSemidefiniteConeScaled <: AbstractSet
    dim::Int
end

dimension(s::Union{PositiveSemidefiniteConeScaled, PositiveSemidefiniteConeTriangle}) = (s.dim * (s.dim + 1)) / 2

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


dimension(s::Union{NonNegative,NonPositive,Zero,SecondOrderCone,Integers,Binaries}) = s.dim

"""
    SOS1(weights::Vector{T}) where T

The set corresponding to the special ordered set (SOS) constraint of type 1. Of the variables in the set, at most one can be nonzero. The ``weights`` induce an ordering of the variables; as such, they should be unique values. The ``k``-th element in the set corresponds to the ``k``-th weight in ``weights``.

See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm) for a description of SOS constraints and their potential uses.
"""
struct SOS1{T} <: AbstractSet
    weights::Vector{T}
end
dimension(s::SOS1) = length(s.weights)

"""
    SOS2(weights::Vector{T}) where T

The set corresponding to the special ordered set (SOS) constraint of type 2. Of the variables in the set, at most two can be nonzero, and if two are nonzero, they must be adjacent in the ordering of the set. The ``weights`` induce an ordering of the variables; as such, they should be unique values. The ``k``-th element in the set corresponds to the ``k``-th weight in ``weights``.

See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm) for a description of SOS constraints and their potential uses.
"""
struct SOS2{T} <: AbstractSet
    weights::Vector{T}
end
dimension(s::SOS2) = length(s.weights)
