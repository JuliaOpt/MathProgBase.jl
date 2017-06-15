

"""
    VariablewiseConstraintReference{T}

A lightweight object used to reference variablewise constraints in a model.
The parameter `T` is the type of set constraint referenced.
"""
struct VariablewiseConstraintReference{T}
    value::UInt64
end

"""
    AffineConstraintReference{T}

A lightweight object used to reference affine-in-set constraints in a model.
The parameter `T` is the type of set constraint referenced.
"""
struct AffineConstraintReference{T}
    value::UInt64
end

"""
    QuadraticConstraintReference{T}

A lightweight object used to reference quadratic-in-set constraints in a model.
The parameter `T` is the type of set constraint referenced.
"""
struct QuadraticConstraintReference{T}
    value::UInt64
end

const ConstraintReference = Union{VariablewiseConstraintReference,AffineConstraintReference,QuadraticConstraintReference}

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
    VariableReference

A lightweight object used to reference variables in a model.
"""
struct VariableReference
    value::UInt64
end

"""
    candelete(m::AbstractMathProgModel, ref::VariableReference)::Bool

Return a `Bool` indicating whether this variable can be removed from the model `m`.
"""
candelete(m::AbstractMathProgModel, ref::VariableReference) = throw(MethodError())

"""
    isvalid(m::AbstractMathProgModel, ref::VariableReference)::Bool

Return a `Bool` indicating whether this reference is valid
for an active variable in the model `m`.
"""
isvalid(m::AbstractMathProgModel, ref::VariableReference) = throw(MethodError())

"""
    delete!(m::AbstractMathProgModel, ref::VariableReference)
    
Delete the referenced variable from the model.

    delete!(m::AbstractMathProgModel, refs::Vector{VariableReference})

Delete the referenced variables in the vector `refs` from the model.
"""
Base.delete!(m::AbstractMathProgModel, ref::VariableReference) = throw(MethodError())
Base.delete!(m::AbstractMathProgModel, refs::Vector{VariableReference}) = throw(MethodError())
