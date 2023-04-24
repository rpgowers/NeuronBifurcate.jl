using Parameters

abstract type ML_Model end
abstract type WB_Model end

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

@with_kw mutable struct MLDS_Param <: ML_Model
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
  ρ::Float64 = 1.0
  λ::Float64 = 100.0
  τδ::Float64 = 10.0
  xin::Float64 = 0.0
end

@with_kw mutable struct MLFDS_Param <: ML_Model
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
  ρ::Float64 = 1.0
  λ::Float64 = 200.0
  τδ::Float64 = 10.0
  L::Float64 = 1000.0
  xin::Float64 = 0.0
end

@with_kw mutable struct MLMDS_Param <: ML_Model
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
  ρ::Float64 = 1.0
  λ::Float64 = 200.0
  τδ::Float64 = 10.0
  L::Float64 = 1000.0
  M::Int64 = 50
  xin::Float64 = 0.0
end

@with_kw mutable struct MLADS_Param <: ML_Model # fully active dendrite and soma
  adims::Int64 = 1 # dimensions of the somatic compartment
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
  ρ::Float64 = 1.0
  λ::Float64 = 200.0
  τδ::Float64 = 10.0
  L::Float64 = 1000.0
  M::Int64 = 50
  xin::Float64 = 0.0
  ϵ::Float64 = 0.1
end

export ML_Model, MLS_Param, MLDS_Param, MLFDS_Param, MLMDS_Param, MLADS_Param

@with_kw mutable struct WBS_Param <: WB_Model
  dims::Int64 = 3 # dimensions of the somatic compartment
  C::Float64 = 1.0
  gL::Float64 = 0.1
  EL::Float64 = -65.0
  gNa::Float64 = 35.0
  ENa::Float64 = 55.0
  gK::Float64 = 9.0
  EK::Float64 = -90.0
  ϕ::Float64 = 5.0
  Iext::Float64 = 0.0
end

@with_kw mutable struct WBDS_Param <: WB_Model
  dims::Int64 = 3 # dimensions of the somatic compartment
  C::Float64 = 1.0
  gL::Float64 = 0.1
  EL::Float64 = -65.0
  gNa::Float64 = 35.0
  ENa::Float64 = 55.0
  gK::Float64 = 9.0
  EK::Float64 = -90.0
  ϕ::Float64 = 5.0
  Iext::Float64 = 0.0
  ρ::Float64 = 1.0
  λ::Float64 = 100.0
  τδ::Float64 = 10.0
  xin::Float64 = 0.0
end

@with_kw mutable struct WBFDS_Param <: WB_Model
  dims::Int64 = 3 # dimensions of the somatic compartment
  C::Float64 = 1.0
  gL::Float64 = 0.1
  EL::Float64 = -65.0
  gNa::Float64 = 35.0
  ENa::Float64 = 55.0
  gK::Float64 = 9.0
  EK::Float64 = -90.0
  ϕ::Float64 = 5.0
  Iext::Float64 = 0.0
  ρ::Float64 = 1.0
  λ::Float64 = 200.0
  τδ::Float64 = 10.0
  L::Float64 = 1000.0
  xin::Float64 = 0.0
end

export WB_Model, WBS_Param, WBDS_Param, WBFDS_Param