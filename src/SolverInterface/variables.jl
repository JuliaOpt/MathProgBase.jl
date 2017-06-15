
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


"""
    addvariables!(m::AbstractMathProgModel, N::Int)

Add `N` scalar variables to the model, returning a vector of variable
references.
"""
function addvariables! end

"""
    addvariable!(m::AbstractMathProgModel, N::Int)

Add a scalar variable to the model, returning a vector of variable
references.
"""
function addvariable! end
