"""
    AbstractSet

Abstract supertype for set objects used to encode constraints.
"""
abstract type AbstractSet end

dimension(s::AbstractSet) = s.dim

"""
    Free(dim)

The set ``\\mathbb{R}^{dim}`` (containing all points) of dimension ``dim``.
"""
struct Zero <: AbstractSet
    dim::Int
end

"""
    Zero(dim)

The set ``\\{ 0 \\}^{dim}`` (containing only the origin) of dimension ``dim``.
"""
struct Zero <: AbstractSet
    dim::Int
end

"""
    NonNegative(dim)

The nonnegative orthant ``\\{ x \\in \\mathbb{R}^{dim} : x \\ge 0 \\}`` of dimension ``dim``.
"""
struct NonNegative <: AbstractSet
    dim::Int
end

"""
    NonPositive(dim)

The nonpositive orthant ``\\{ x \\in \\mathbb{R}^{dim} : x \\le 0 \\}`` of dimension ``dim``.
"""
struct NonPositive <: AbstractSet
    dim::Int
end

"""
    Interval(lower,upper)

The box ``[lower, upper] \\subseteq \\mathbb{R}^{dim}`` where ``lower`` and ``upper`` are vectors of dimension ``dim``. If ``lower`` or ``upper`` is all ``-Inf`` or ``Inf``, the set is interpreted as a one-sided interval.
"""
struct Interval{T <: Real} <: AbstractSet
    lower::T
    upper::T
end

dimension(s::Interval) = length(s.lower)

"""
    SecondOrderCone(dim)

The second-order cone (or Lorenz cone) ``\\{ (t,x) \\in \\mathbb{R}^{dim} : t \\ge || x ||_2 \\}`` of dimension ``dim``.
"""
struct SecondOrderCone <: AbstractSet
    dim::Int
end

"""
    ExponentialCone()

The 3-dimensional exponential cone ``\\{ (x,y,z) \\in \\mathbb{R}^3 : y \\exp{x/y} \\le z, y > 0 \\}``.
"""
struct ExponentialCone <: AbstractSet
end

"""
    DualExponentialCone()

The 3-dimensional dual exponential cone ``\\{ (u,v,w) \\in \\mathbb{R}^3 : -u \\exp{v/u} \\le \\exp{1} w, u < 0 \\}``.
"""
struct DualExponentialCone <: AbstractSet
end

"""
    PowerCone(a)

The 3-dimensional power cone ``\\{ (x,y,z) \\in \\mathbb{R}^3 : x^{a} y^{1-a} >= |z|, x \\ge 0, y \\ge 0 \\}`` with parameter ``a``.
"""
struct PowerCone{T <: Real} <: AbstractSet
    a::T
end

"""
    DualPowerCone(a)

The 3-dimensional power cone ``\\{ (u,v,w) \\in \\mathbb{R}^3 : (u/a)^a (v/(1-a))^{1-a} >= |w|, u \\ge 0, v \\ge 0 \\}`` with parameter ``a``.
"""
struct DualPowerCone{T <: Real} <: AbstractSet
    a::T
end

dimension(s::Union{ExponentialCone, DualExponentialCone, PowerCone, DualPowerCone}) = 3

"""
    PositiveSemidefiniteConeTriangle(dim)

The cone of symmetric ``dim \\times dim`` matrices that are positive semidefinite. The dimension of the cone is ``dim (dim + 1)/2`` since the matrices are symmetric. The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column). The scalar product is the sum of the pairwise product of the diagonal entries plus twice the sum of the pairwise product of the upper diagonal entries.

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
    sidedim::Int
end

"""
    PositiveSemidefiniteConeScaled(n)

The cone of symmetric ``dim \\times dim`` matrices that are positive semidefinite. The dimension of the cone is ``dim (dim + 1)/2`` since the matrices are symmetric. The entries of the upper triangular part of the matrix are given row by row (or equivalently, the entries of the lower triangular part are given column by column). The off-diagonal entries of the matrices of both the cone and its dual are scaled by ``\\sqrt{2}`` and the scalar product is simply the sum of the pairwise product of the entries.

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
    sidedim::Int
end

"""
    Integers(dim)

The set of integers ``\\mathbb{Z}^{dim}``.
"""
struct Integers <: AbstractSet
    dim::Int
end

"""
    Binaries(dim)

The set of binary vectors ``\\{ 0, 1 \\}^{dim}``.
"""
struct Binaries <: AbstractSet
    dim::Int
end

"""
    SOS1(weights)

The set corresponding to the special ordered set (SOS) constraint of type 1. Of the variables in the set, at most one can be nonzero. The ``weights`` induce an ordering of the variables; as such, they should be unique values. The ``k``-th element in the set corresponds to the ``k``-th weight in ``weights``. See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm) for a description of SOS constraints and their potential uses.
"""
struct SOS1{T <: Real} <: AbstractSet
    weights::Vector{T}
end

"""
    SOS2(weights)

The set corresponding to the special ordered set (SOS) constraint of type 2. Of the variables in the set, at most two can be nonzero, and if two are nonzero, they must be adjacent in the ordering of the set. The ``weights`` induce an ordering of the variables; as such, they should be unique values. The ``k``-th element in the set corresponds to the ``k``-th weight in ``weights``. See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm) for a description of SOS constraints and their potential uses.
"""
struct SOS2{T <: Real} <: AbstractSet
    weights::Vector{T}
end

dimension(s::Union{SOS1, SOS2}) = length(s.weights)
