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

function Ia((n,v),args)
  @unpack gs, Es, Am, Δm, gf, Ef = args
  gs*n*(Es-v)+gf*a∞(v,Am,Δm)*(Ef-v)
end

function X∞(v,args::ML_Model)
  [a∞(v,args.An,args.Δn), v]
end

function dA∞(v,args::ML_Model)
  da∞(v,args.An,args.Δn)
end

function τa(v,args::ML_Model)
  1/ψa(v,args.ϕ,args.An,args.Δn)
end

function F(x, args::MLS_Param) # note that this does not currently include Iext in soma_voltage
  [ML_ncurrent(x,args), soma_voltage(x,args)]
end

# function I∞(v,args::Union{MLS_Param})
#   return soma_voltage(X∞(v,args),args)*args.C+args.Iext
#   # return soma_voltage((a∞(v,args.An,args.Δn),v),args)*args.C+args.Iext
# end

# function I∞(v, args::Union{MLDS_Param})
#   @unpack C, ρ, xin, λ, gL, EL, Iext = args
#   return soma_voltage(X∞(v,args),args)*C+ρ*gL*(EL-v)+Iext*exp(-xin/λ)
#   # return soma_voltage((a∞(v,An,Δn),v),args)*C+ρ*gL*(EL-v)+Iext*exp(-xin/λ)
# end

function ∂Ia∂n(v,args::ML_Model)
  @unpack gs, Es = args
  gs*(Es-v)
end

function ∂Ia∂v(v,args::ML_Model)
  @unpack gs, An, Δn, gf, Ef, Am, Δm = args
  -gs*a∞(v,An,Δn)-gf*a∞(v,Am,Δm)+gf*da∞(v,Am,Δm)*(Ef-v)
end

# function bt(args::Union{MLS_Param})
#   @unpack C, gs, Es, ϕ, An, Δn = args
#   f(v) = 1+∂Ia∂n(v,args)*da∞(v,An,Δn)/(C*ψa(v,ϕ,An,Δn))
#   vbt = find_zeros(v->f(v), -80.0, 20.0)
#   args_temp = @set args.gL = 0.0
#   g(v) = I∞(v,args_temp)
#   gbt = ForwardDiff.derivative.(g,vbt)
#   Ibt = zeros(length(gbt))
#   for i in eachindex(vbt)
#     args_temp = setproperties(args, (gL=gbt[i], Iext=0.0))
#     Ibt[i] = -I∞(vbt[i],args_temp)/Iscale(args_temp)
#   end
#   return vbt, Ibt, gbt
# end

# function bt(args::Union{MLDS_Param})
#   @unpack C, gL, An, Δn, ϕ, τδ = args
#   f(v) = C+0.5*τδ*(-gL+∂Ia∂v(v,args))+(0.5*τδ+1/ψa(v,ϕ,An,Δn))*∂Ia∂n(v,args)*da∞(v,An,Δn)
#   ρ(v) = (-gL+∂Ia∂v(v,args)+∂Ia∂n(v,args)*da∞(v,An,Δn))/gL
#   vbt = find_zeros(v -> f(v), -80.0, 40.0)
#   ρbt = ρ.(vbt)
#   Ibt = zeros(length(vbt))
#   for i in eachindex(Ibt)
#     args_temp = setproperties(args, (ρ = ρbt[i], Iext=0.0))
#     Ibt[i] = -I∞(vbt[i],args_temp)/Iscale(args_temp)
#   end
#   return vbt, Ibt, ρbt
# end

function hopf(args::MLS_Param)
  @unpack An, Δn, ϕ, gL, C = args
  f(v) = -gL+∂Ia∂v(v,args)-C*ψa(v,ϕ,An,Δn)
  vh = find_zeros(v->f(v), -80.0, 30.0)
  ω2(v) = -ψa(v,ϕ,An,Δn)*(ψa(v,ϕ,An,Δn)+∂Ia∂n(v,args)*da∞(v,An,Δn)/C)
  args_temp = @set args.Iext = 0.0
  Ih = zeros(length(vh))
  ωh = zeros(length(vh))
  for i in eachindex(vh)
    if ω2(vh[i]) < 0
      Ih[i] = NaN
      ωh[i] = NaN
    else
      Ih[i] = -I∞(vh[i],args_temp)
      ωh[i] = sqrt(ω2(vh[i]))
    end
  end
  return vh, Ih, ωh
end