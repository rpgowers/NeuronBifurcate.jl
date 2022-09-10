using Parameters

abstract type ML_Model end

@with_kw mutable struct MLS_Param <: ML_Model
    dims::Int64 = 2 # dimensions of the somatic compartment
    C::Float64 = 20.0
    Iext::Float64 = 0.0
    gL::Float64 = 2.0
    EL::Float64 = -60.0
    gf::Float64 = 4.0
    Ef::Float64 = 120.0
    gs::Float64 = 8.0
    Es::Float64 = -80.0
    Am::Float64 = -1.2
    Δm::Float64 = 9.0
    An::Float64 = 12.0
    Δn::Float64 = 8.7
    ϕ::Float64 = 1/15
  end
  export ML_Model, MLS_Param