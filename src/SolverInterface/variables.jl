"""
    addvariables!(m::AbstractMathProgModel, N::Int)

Add `N` scalar variables to the model, returning a vector of variable
references.
"""
function addvariables! end

"""
    addvariable!(m::AbstractMathProgModel, N::Int)

Add a scalar variable to the model, returning variable reference.
"""
function addvariable! end
