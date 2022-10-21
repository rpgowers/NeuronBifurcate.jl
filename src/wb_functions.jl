αm(v) = (v == -35 ? 1.0 : (-(v+35)/10)/(exp(-(v+35)/10)-1))
βm(v) = 4*exp(-(v+60)/18)
minf(v) = αm(v)/(αm(v)+βm(v))

αn(v) = (v == -34 ? 0.1 : (-(v+34)/100)/(exp(-(v+34)/10)-1) )
βn(v) = exp(-(v+44)/80)/8
dαn(v) = (v==-34 ? 1/200 : -0.01*( exp(-(v+34)/10)*(1+(v+34)/10)-1 )/((exp(-(v+34)/10)-1)^2) )
dβn(v) = -exp(-(v+44)/80)/640

τn(v, ϕ) = 1/(ϕ*(αn(v)+βn(v)))
ninf(v) = αn(v)/(αn(v)+βn(v))
dninf(v) = (dαn(v)*βn(v)-αn(v)*dβn(v))/(αn(v)+βn(v))^2

αh(v) = 0.07*exp(-(v+58)/20)
βh(v) = 1/(1+exp(-(v+28)/10))
dαh(v) = -0.07/20 * exp(-(v+58)/20)
dβh(v) = 0.1*exp(-(v+28)/10)/(1+exp(-(v+28)/10))^2

τh(v, ϕ) = 1/(ϕ*(αh(v)+βh(v)))
hinf(v) = αh(v)/(αh(v)+βh(v))
dhinf(v) = (dαh(v)*βh(v)-αh(v)*dβh(v))/(αh(v)+βh(v))^2

function WB_neq((n,h,v),args)
  @unpack ϕ = args
  return (ninf(v)-n)/τn(v, ϕ)
end
  
function WB_heq((n,h,v),args)
  @unpack ϕ = args
  return (hinf(v)-h)/τh(v, ϕ)
end
  
function soma_voltage((n,h,v),args::WBS_Param)
  @unpack C, gL, EL, gNa, ENa, gK, EK, Iext = args
  return (gL*(EL-v)+gNa*minf(v)^3*h*(ENa-v)+gK*n^4*(EK-v)+Iext)/C
end

function F(x, args::WBS_Param)
  [WB_neq(x,args), WB_heq(x,args), soma_voltage(x,args)]
end

function I∞(v,args::WBS_Param)
  return soma_voltage((ninf(v),hinf(v),v),args)*args.C
end

function ∂Ia∂n(v,args::WB_Model) # as evaluated at equilibrium
  @unpack C, EK, gK = args
  4*gK*ninf(v)^3*(EK-v)
end

function ∂Ia∂h(v,args::WB_Model) # as evaluated at equilibrium
  @unpack C, ENa, gNa = args
  gNa*minf(v)^3*(ENa-v)
end

function bt(args::WBS_Param)
  @unpack C, gNa, ENa, gK, EK, ϕ = args
  vbt = find_zeros(v-> 1+τn(v,ϕ)*∂Ia∂n(v,args)*dninf(v)/C+τh(v,ϕ)*∂Ia∂h(v,args)*dhinf(v)/C , -80.0, 40.0)
  args_temp = @set args.gL = 0.0
  f2(v) = I∞(v,args_temp)
  gbt = ForwardDiff.derivative.(f2, vbt)
  Ibt = zeros(length(gbt))
  for i in eachindex(vbt)
    args_temp = setproperties(args, (gL=gbt[i], Iext=0.0))
    Ibt[i] = -I∞(vbt[i],args_temp)/Iscale(args_temp)
  end
  return vbt, Ibt, gbt
end