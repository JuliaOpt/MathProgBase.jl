# Methods for the Conic interface

abstract type AbstractConicModel <: AbstractMathProgModel end
export AbstractConicModel

@define_interface begin
    ConicModel
    getdual
    getvardual
    supportedcones
    setbvec!
end
