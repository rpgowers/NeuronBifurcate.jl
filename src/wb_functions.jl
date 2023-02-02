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

function Ia((n,h,v), args::WBS_Param)
  @unpack gNa, ENa, gK, EK = args
  gNa*minf(v)^3*h*(ENa-v)+gK*n^4*(EK-v)
end

function X∞(v,args::WBS_Param)
  [ninf(v), hinf(v), v]
end

function dA∞(v,args::WBS_Param)
  [dninf(v), dhinf(v)]
end

function τa(v,args::WBS_Param)
  [τn(v, args.ϕ), τh(v, args.ϕ)]
end
  
function F(x, args::WBS_Param)
  [WB_neq(x,args), WB_heq(x,args), soma_voltage(x,args)]
end

function ∂Ia∂n(v,args::WB_Model) # as evaluated at equilibrium
  @unpack C, EK, gK = args
  4*gK*ninf(v)^3*(EK-v)
end

function ∂Ia∂h(v,args::WB_Model) # as evaluated at equilibrium
  @unpack C, ENa, gNa = args
  gNa*minf(v)^3*(ENa-v)
end