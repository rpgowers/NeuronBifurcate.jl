αm(v) = vtrap(-0.1, v+35, -10)
# αm(v) = (v == -35 ? 1.0 : (-(v+35)/10)/(exp(-(v+35)/10)-1))
βm(v) = 4*exp(-(v+60)/18)
dαm(v) = (v == -35 ? 0.05 : -0.1*( exp(-(v+35)/10)*(1+(v+35)/10)-1 )/((exp(-(v+35)/10)-1)^2) )
dβm(v) = -4/18 * exp(-(v+60)/18)
d2αm(v) = (v==-35 ? 1/600 : -exp(-(v+35)/10)*(exp(-(v+35)/10)*(2+(v+35)/10)+(v+35)/10-2)/(100*(exp(-(v+35)/10)-1)^3))
d2βm(v) = exp(-(v+60)/18)/81
d3αm(v) = (v == - 35 ? 0.0 : -exp(-(v+35)/10)*(4*(v+35)/10*exp(-(v+35)/10)+((v+35)/10+3)*exp(-2*(v+35)/10)-3+(v+35)/10 )/(1000*(exp(-(v+35)/10)-1)^4) )
d3βm(v) = -exp(-(v+60)/18)/1458

minf(v) = αm(v)/(αm(v)+βm(v))
dminf(v) = (dαm(v)*βm(v)-αm(v)*dβm(v))/(αm(v)+βm(v))^2
d2minf(v) =( (αm(v)+βm(v))*(d2αm(v)*βm(v)-αm(v)*d2βm(v))-2(dαm(v)*βm(v)-αm(v)*dβm(v))*(dαm(v)+dβm(v)) )/(αm(v)+βm(v))^3

function d3minf(v)
  t1 = (αm(v)+βm(v))*(d3αm(v)*βm(v)+d2αm(v)*dβm(v)-αm(v)*d3βm(v)-dαm(v)*d2βm(v))
  t2 = -(d2αm(v)*βm(v)-αm(v)*d2βm(v))*(dαm(v)+dβm(v))-2*(dαm(v)*βm(v)-αm(v)*dβm(v))*(d2αm(v)+d2βm(v))
  t3 = -3*(dαm(v)+dβm(v))*( (αm(v)+βm(v))*(d2αm(v)*βm(v)-αm(v)*d2βm(v))-2(dαm(v)*βm(v)-αm(v)*dβm(v))*(dαm(v)+dβm(v)) )
  return ((αm(v)+βm(v))*(t1+t2)+t3)/(αm(v)+βm(v))^4
end

αh(v) = 0.07*exp(-(v+58)/20)
βh(v) = 1/(1+exp(-(v+28)/10))
dαh(v) = -0.07/20 * exp(-(v+58)/20)
dβh(v) = 0.1*exp(-(v+28)/10)/(1+exp(-(v+28)/10))^2
d2αh(v) = 0.07/400 * exp(-(v+58)/20)
d2βh(v) = ( exp(-2*(v+28)/10)-exp(-(v+28)/10) )/(100*(1+exp(-(v+28)/10))^3)
d3αh(v) = -0.07/8000 * exp(-(v+58)/20)
d3βh(v) = (exp(-3*(v+28)/10)-4*exp(-2*(v+28)/10)+exp(-(v+28)/10))/(1000*(exp(-(v+28)/10)+1)^4)

τh(v, ϕ) = 1/(ϕ*(αh(v)+βh(v)))
dτh(v, ϕ) = 1/(ϕ*(dαh(v)+dβh(v))) # this is 1/d(1/τh)/dv
d2τh(v, ϕ) = 1/(ϕ*(d2αh(v)+d2βh(v))) # this is 1/d2(1/τh)/dv2

hinf(v) = αh(v)/(αh(v)+βh(v))
dhinf(v) = (dαh(v)*βh(v)-αh(v)*dβh(v))/(αh(v)+βh(v))^2
d2hinf(v) =( (αh(v)+βh(v))*(d2αh(v)*βh(v)-αh(v)*d2βh(v))-2(dαh(v)*βh(v)-αh(v)*dβh(v))*(dαh(v)+dβh(v)) )/(αh(v)+βh(v))^3

function d3hinf(v)
  t1 = (αh(v)+βh(v))*(d3αh(v)*βh(v)+d2αh(v)*dβh(v)-αh(v)*d3βh(v)-dαh(v)*d2βh(v))
  t2 = -(d2αh(v)*βh(v)-αh(v)*d2βh(v))*(dαh(v)+dβh(v))-2*(dαh(v)*βh(v)-αh(v)*dβh(v))*(d2αh(v)+d2βh(v))
  t3 = -3*(dαh(v)+dβh(v))*( (αh(v)+βh(v))*(d2αh(v)*βh(v)-αh(v)*d2βh(v))-2(dαh(v)*βh(v)-αh(v)*dβh(v))*(dαh(v)+dβh(v)) )
  return ((αh(v)+βh(v))*(t1+t2)+t3)/(αh(v)+βh(v))^4
end

# αn(v) = (v == -34 ? 0.1 : (-(v+34)/100)/(exp(-(v+34)/10)-1) )
αn(v) = vtrap(-0.01, v+34, -10)

βn(v) = exp(-(v+44)/80)/8
dαn(v) = (v==-34 ? 1/200 : -0.01*( exp(-(v+34)/10)*(1+(v+34)/10)-1 )/((exp(-(v+34)/10)-1)^2) )
dβn(v) = -exp(-(v+44)/80)/640
d2αn(v) = (v==-34 ? 1/6000 : -exp(-(v+34)/10)*(exp(-(v+34)/10)*(2+(v+34)/10)+(v+34)/10-2)/(1000*(exp(-(v+34)/10)-1)^3))
d2βn(v) = exp(-(v+44)/80)/51200
d3αn(v) = (v==-34 ? 0.0 : - (exp(-(v+34)/10)*((v+34)/10-3)+4*(v+34)/10*exp(-(2*(v+34)/10))+exp(-(3*(v+34)/10))*((v+34)/10+3))/(10000*(exp(-(v+34)/10)-1)^4) )
d3βn(v) = -exp(-(v+44)/80)/4096000

τn(v, ϕ) = 1/(ϕ*(αn(v)+βn(v)))
dτn(v, ϕ) = 1/(ϕ*(dαn(v)+dβn(v))) # this is 1/d(1/τn)/dv
d2τn(v, ϕ) = 1/(ϕ*(d2αn(v)+d2βn(v))) # this is 1/d2(1/τn)/dv2

ninf(v) = αn(v)/(αn(v)+βn(v))
dninf(v) = (dαn(v)*βn(v)-αn(v)*dβn(v))/(αn(v)+βn(v))^2
d2ninf(v) = ( (αn(v)+βn(v))*(d2αn(v)*βn(v)-αn(v)*d2βn(v))-2(dαn(v)*βn(v)-αn(v)*dβn(v))*(dαn(v)+dβn(v)) )/(αn(v)+βn(v))^3

function d3ninf(v)
  t1 = (αn(v)+βn(v))*(d3αn(v)*βn(v)+d2αn(v)*dβn(v)-αn(v)*d3βn(v)-dαn(v)*d2βn(v))
  t2 = -(d2αn(v)*βn(v)-αn(v)*d2βn(v))*(dαn(v)+dβn(v))-2*(dαn(v)*βn(v)-αn(v)*dβn(v))*(d2αn(v)+d2βn(v))
  t3 = -3*(dαn(v)+dβn(v))*( (αn(v)+βn(v))*(d2αn(v)*βn(v)-αn(v)*d2βn(v))-2(dαn(v)*βn(v)-αn(v)*dβn(v))*(dαn(v)+dβn(v)) )
  return ((αn(v)+βn(v))*(t1+t2)+t3)/(αn(v)+βn(v))^4
end

export ninf
export hinf

## should be able to strip out most of the above eventually

function WB_neq((n,h,v),args)
  @unpack ϕ = args
  return (ninf(v)-n)/τn(v, ϕ)
end
  
function WB_heq((n,h,v),args)
  @unpack ϕ = args
  return (hinf(v)-h)/τh(v, ϕ)
end

function Ia((n,h,v), args::WB_Model)
  @unpack gNa, ENa, gK, EK = args
  gNa*minf(v)^3*h*(ENa-v)+gK*n^4*(EK-v)
end

function X∞(v,args::WB_Model)
  [ninf(v), hinf(v), v]
end

function dA∞(v,args::WB_Model)
  [dninf(v), dhinf(v)]
end

function τa(v,args::WB_Model)
  [τn(v, args.ϕ), τh(v, args.ϕ)]
end
  
function F(x, args::WB_Model)
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

function ∂f∂v(v,args::WB_Model)
  @unpack C, gL, gNa, ENa, gK = args
  return (-gL-gK*ninf(v)^4+gNa*hinf(v)*(-minf(v)^3+3*minf(v)^2*dminf(v)*(ENa-v)))/C
end

function ∂f∂n(v,args::WB_Model)
  @unpack C, EK, gK = args
  4*gK*ninf(v)^3*(EK-v)/C
end

function ∂f∂h(v,args::WB_Model)
  @unpack C, ENa, gNa = args
  gNa*minf(v)^3*(ENa-v)/C
end

function hopf(args::WBS_Param)
  @unpack C, ϕ = args

  a(v) = τn(v,ϕ)^2*τh(v,ϕ)^2
  b(v) = τn(v,ϕ)*τh(v,ϕ)^2*∂f∂n(v,args)*dninf(v)+τh(v,ϕ)*τn(v,ϕ)^2*∂f∂h(v,args)*dhinf(v)+τn(v,ϕ)^2+τh(v,ϕ)^2
  c(v) = τn(v,ϕ)*∂f∂n(v,args)*dninf(v)+τh(v,ϕ)*∂f∂h(v,args)*dhinf(v)+1
  Ω(v) = (-b(v)+sqrt(b(v)^2-4*a(v)*c(v)))/(2*a(v))
  F(v) = (1+Ω(v)*τh(v,ϕ)^2)*∂f∂n(v,args)*dninf(v)+(1+Ω(v)*τn(v,ϕ)^2)*∂f∂h(v,args)*dhinf(v)+(1+Ω(v)*τn(v,ϕ)^2)*(1+Ω(v)*τh(v,ϕ)^2)*∂f∂v(v,args)
  vin = find_zeros(v->F(v), -80.0, 60.0) # initial values of vh
  vh = NaN.*ones(length(vin))
  ωh = NaN.*ones(length(vin))
  Ih = NaN.*ones(length(vin))
  for i in eachindex(vin)
    Ωin = Ω(vin[i])
    args_temp = @set args.Iext = 0.0
    if Ωin > 0.0
      vh[i] = vin[i]
      ωh[i] = sqrt(Ωin)
      Ih[i] = -I∞(vh[i], args_temp)
    end
  end
  return vh, Ih, ωh
end

function Bexp_old(x,y,v,args::WB_Model)
  @unpack C, gNa, ENa, gK, EK, ϕ = args
  n = ninf(v)
  h = hinf(v)

  dFdv2 = gNa*h*(-6*minf(v)^2*dminf(v)+6*minf(v)*dminf(v)^2*(ENa-v)+3*minf(v)^2*d2minf(v)*(ENa-v))/C
  dFdvdn = -4*gK*n^3/C
  dFdn2 = 12*gK*n^2*(EK-v)/C
  dFdvdh = gNa*(-minf(v)^3+3*minf(v)^2*dminf(v)*(ENa-v) )/C
  Bv = dFdv2*x[3]*y[3]+dFdvdn*(x[3]*y[1]+x[1]*y[3])+dFdvdh*(x[3]*y[2]+x[2]*y[3])+dFdn2*x[1]*y[1]

  dNdv2 = d2ninf(v)/τn(v,ϕ)+2*dninf(v)/dτn(v,ϕ)
  dNdvdn = -1/dτn(v,ϕ)
  Bn = dNdv2*x[3]*y[3]+dNdvdn*(x[3]*y[1]+x[1]*y[3])

  dHdv2 = d2hinf(v)/τh(v,ϕ)+2*dhinf(v)/dτh(v,ϕ)
  dHdvdn = -1/dτh(v,ϕ)
  Bh = dHdv2*x[3]*y[3]+dHdvdn*(x[3]*y[2]+x[2]*y[3])
  return [Bn, Bh, Bv]
end

function Cexp_old(x,y,z,v,args::WB_Model)
  @unpack C, gNa, ENa, gK, EK, ϕ = args
  n = ninf(v)
  h = hinf(v)

  dFdv3 = gNa*h*(-18*minf(v)*dminf(v)^2-9*minf(v)^2*d2minf(v)+6*dminf(v)^3*(ENa-v)+18*minf(v)*dminf(v)*d2minf(v)*(ENa-v)+3*minf(v)^2*d3minf(v)*(ENa-v))/C
  dFdvdn2 = -12*gK*n^2/C
  dFdn3 = 24*gK*n*(EK-v)/C
  dFdv2dh = gNa*(-6*minf(v)^2*dminf(v)+6*minf(v)*dminf(v)^2*(ENa-v)+3*minf(v)^2*d2minf(v)*(ENa-v))/C

  dNdv2dn = -1/d2τn(v,ϕ)
  dNdv3 = d3ninf(v)/τn(v,ϕ)+3*d2ninf(v)/dτn(v,ϕ)+3*dninf(v)/d2τn(v,ϕ)

  dHdv2dh = -1/d2τh(v,ϕ)
  dHdv3 = d3hinf(v)/τh(v,ϕ)+3*d2hinf(v)/dτh(v,ϕ)+3*dhinf(v)/d2τh(v,ϕ)

  Cv = dFdv3*x[3]*y[3]*z[3]+dFdvdn2*(x[3]*y[1]*z[1]+x[1]*y[3]*z[1]+x[1]*y[1]*z[3])+dFdn3*x[1]*y[1]*z[1]+dFdv2dh*(x[3]*y[3]*z[2]+x[3]*y[2]*z[3]+x[2]*y[3]*z[3])
  Cn = dNdv2dn*(x[3]*y[3]*z[1]+x[3]*y[1]*z[3]+x[1]*y[3]*z[3])+dNdv3*x[3]*y[3]*z[3]
  Ch = dHdv2dh*(x[3]*y[3]*z[2]+x[3]*y[2]*z[3]+x[2]*y[3]*z[3])+dHdv3*x[3]*y[3]*z[3]
  return [Cn, Ch, Cv]
end