# Methods for the Conic interface

@compat abstract type AbstractConicModel <: AbstractMathProgModel end
export AbstractConicModel

@define_interface begin
    ConicModel
    getdual
    getvardual
    supportedcones
    setbvec!
end
