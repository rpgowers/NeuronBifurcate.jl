function a∞(v,A,Δ)
  return 1/( 1+exp( -(v-A)/Δ ) )
end
export a∞

function da∞(v,A,Δ)
  return exp( -(v-A)/Δ )/Δ/( 1+exp( -(v-A)/Δ ) )^2
end

function ψa(v,ϕ,A,Δ)
  return ϕ*cosh( (v-A)/(4*Δ) )
end
export ψa

function ML_ncurrent((n,v),args)
  @unpack An, Δn, ϕ = args
  return (a∞(v,An,Δn)-n)*ψa(v,ϕ,An,Δn)
end
export ML_ncurrent

function soma_voltage((n,v),args::MLS_Param)
  @unpack C, EL, gL, Ef, gf, Am, Δm, Es, gs, Iext = args
  return (gL*(EL-v)+gs*n*(Es-v)+gf*a∞(v,Am,Δm)*(Ef-v)+Iext)/C
end

function F(x, args::MLS_Param)
  [ML_ncurrent(x,args), soma_voltage(x,args)]
end

function I∞(v,args::Union{MLS_Param})
  return soma_voltage((a∞(v,args.An,args.Δn),v),args)*args.C
end

function ∂Ia∂n(v,args::ML_Model)
  @unpack gs, Es = args
  gs*(Es-v)
end

function bt(args::Union{MLS_Param})
  @unpack C, gs, Es, ϕ, An, Δn = args
  f(v) = 1+∂Ia∂n(v,args)*da∞(v,An,Δn)/(C*ψa(v,ϕ,An,Δn))
  vbt = find_zeros(v->f(v), -80.0, 20.0)
  args_temp = @set args.gL = 0.0
  g(v) = I∞(v,args_temp)
  gbt = ForwardDiff.derivative.(g,vbt)
  Ibt = zeros(length(gbt))
  for i in eachindex(vbt)
    args_temp = setproperties(args, (gL=gbt[i], Iext=0.0))
    Ibt[i] = -I∞(vbt[i],args_temp)/Iscale(args_temp)
  end
  return vbt, Ibt, gbt
end