"""
    addvariables!(m::AbstractMathProgModel, N::Int)::Vector{VariableReference}

Add `N` scalar variables to the model, returning a vector of variable
references.
"""
function addvariables! end

"""
    addvariable!(m::AbstractMathProgModel)::VariableReference

Add a scalar variable to the model, returning a variable reference.
"""
function addvariable! end
