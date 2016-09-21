# Methods for the Conic interface

abstract AbstractConicModel <: AbstractMathProgModel
export AbstractConicModel

@define_interface begin
    ConicModel
    getdual
    getvardual
    supportedcones
    setbvec!
end
