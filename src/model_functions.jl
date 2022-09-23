include("model_arguments.jl")

using Roots

function a∞(v,A,Δ)
  return 1/( 1+exp( -(v-A)/Δ ) )
end
export a∞

function ψa(v,ϕ,A,Δ)
  return ϕ*cosh( (v-A)/(4*Δ) )
end
export ψa

function ML_ncurrent((n,v),args)
  @unpack An, Δn, ϕ = args
  return (a∞(v,An,Δn)-n)*ψa(v,ϕ,An,Δn)
end
export ML_ncurrent

function MLS_voltage((n,v),args)
  @unpack C, EL, gL, Ef, gf, Am, Δm, Es, gs, Iext = args
  return (gL*(EL-v)+gs*n*(Es-v)+gf*a∞(v,Am,Δm)*(Ef-v)+Iext)/C
end
export MLS_voltage

function F(x, args::MLS_Param)
  [ML_ncurrent(x,args), MLS_voltage(x,args)]
end
export F

function I∞(v,args::Union{MLS_Param})
  return MLS_voltage((a∞(v,args.An,args.Δn),v),args)*args.C
end
export I∞

function vfps(args)
  vfp = find_zeros(v->I∞(v,args), -120.0, 50.0)
end
export vfps