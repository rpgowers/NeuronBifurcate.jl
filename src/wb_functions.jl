αm(v) = vtrap(-0.1, v+35, -10)
βm(v) = 4*exp(-(v+60)/18)
dαm(v) = (v == -35 ? 0.05 : -0.1*( exp(-(v+35)/10)*(1+(v+35)/10)-1 )/((exp(-(v+35)/10)-1)^2) )
dβm(v) = -4/18 * exp(-(v+60)/18)

minf(v) = αm(v)/(αm(v)+βm(v))
dminf(v) = (dαm(v)*βm(v)-αm(v)*dβm(v))/(αm(v)+βm(v))^2

αh(v) = 0.07*exp(-(v+58)/20)
βh(v) = 1/(1+exp(-(v+28)/10))
dαh(v) = -0.07/20 * exp(-(v+58)/20)
dβh(v) = 0.1*exp(-(v+28)/10)/(1+exp(-(v+28)/10))^2

τh(v, ϕ) = 1/(ϕ*(αh(v)+βh(v)))

hinf(v) = αh(v)/(αh(v)+βh(v))
dhinf(v) = (dαh(v)*βh(v)-αh(v)*dβh(v))/(αh(v)+βh(v))^2

αn(v) = vtrap(-0.01, v+34, -10)
βn(v) = exp(-(v+44)/80)/8
dαn(v) = (v==-34 ? 1/200 : -0.01*( exp(-(v+34)/10)*(1+(v+34)/10)-1 )/((exp(-(v+34)/10)-1)^2) )
dβn(v) = -exp(-(v+44)/80)/640

τn(v, ϕ) = 1/(ϕ*(αn(v)+βn(v)))

ninf(v) = αn(v)/(αn(v)+βn(v))
dninf(v) = (dαn(v)*βn(v)-αn(v)*dβn(v))/(αn(v)+βn(v))^2

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
  
function Fσ_dyn(x, args::WB_Model)
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

# function hopf_old(args::WBDS_Param; v0 =-10.0, ω0=0.0025)
#   @unpack C, ϕ, τδ, ρ, gL = args
#   function WBC1_hopftest!(F, (v, ω))
#     F[1] = ∂f∂n(v,args)*dninf(v)/(1+ω^2*τn(v,ϕ)^2)+∂f∂h(v,args)*dhinf(v)/(1+ω^2*τh(v,ϕ)^2)+∂f∂v(v,args)-ρ*gL*z(ω,τδ)/(2*C)
#     F[2] = τn(v,ϕ)*∂f∂n(v,args)*dninf(v)/(1+ω^2*τn(v,ϕ)^2)+τh(v,ϕ)*∂f∂h(v,args)*dhinf(v)/(1+ω^2*τh(v,ϕ)^2)+1+ρ*gL*u(ω,τδ)/(2*C*ω)
#   end
#   output = nlsolve(WBC1_hopftest!, [v0;ω0], autodiff = :forward, ftol=1e-8)
#   vh = output.zero[1]
#   ωh = output.zero[2]
#   args_temp = @set args.Iext = 0.0
#   Ih = -I∞(vh, args_temp)
#   return vh, Ih, ωh
# end
# export hopf_old

function l1_coeff_general(p,q,B,C,ω,v,J,args)
  @unpack dims = args
  g21 = p'C(q,q,conj.(q),v,args)-2*p'B(q,inv(J)*B(q,conj.(q),v,args),v,args)
  rest = p'B(conj.(q),inv(2im*ω*I(dims)-J)*B(q,q,v,args),v,args)
  return real(g21+rest)/(2*ω)
end

function hopf_stability(v, ω, args::WBS_Param)
  J = Jacobian(X∞(v,args),args) # check this gives the expected eigenvectors and eigenvalues
  q = eigvecs(Array(J))[:,end]
  p = eigvecs(Array(J'))[:,end-1]
  k = p'q
  if ω == 0 # if no frequency is given
    ω = imag(eigvals(Array(J))[end])
  end
  pn = p./conj(k) # this is to normalise the inner product to
  return l1_coeff_general(pn,q,Bexp,Cexp,ω,v,J,args)
end

# function Bexp_old(x,y,v,args::WB_Model)
#   @unpack C, gNa, ENa, gK, EK, ϕ = args
#   n = ninf(v)
#   h = hinf(v)

#   dFdv2 = gNa*h*(-6*minf(v)^2*dminf(v)+6*minf(v)*dminf(v)^2*(ENa-v)+3*minf(v)^2*d2minf(v)*(ENa-v))/C
#   dFdvdn = -4*gK*n^3/C
#   dFdn2 = 12*gK*n^2*(EK-v)/C
#   dFdvdh = gNa*(-minf(v)^3+3*minf(v)^2*dminf(v)*(ENa-v) )/C
#   Bv = dFdv2*x[3]*y[3]+dFdvdn*(x[3]*y[1]+x[1]*y[3])+dFdvdh*(x[3]*y[2]+x[2]*y[3])+dFdn2*x[1]*y[1]

#   dNdv2 = d2ninf(v)/τn(v,ϕ)+2*dninf(v)/dτn(v,ϕ)
#   dNdvdn = -1/dτn(v,ϕ)
#   Bn = dNdv2*x[3]*y[3]+dNdvdn*(x[3]*y[1]+x[1]*y[3])

#   dHdv2 = d2hinf(v)/τh(v,ϕ)+2*dhinf(v)/dτh(v,ϕ)
#   dHdvdn = -1/dτh(v,ϕ)
#   Bh = dHdv2*x[3]*y[3]+dHdvdn*(x[3]*y[2]+x[2]*y[3])
#   return [Bn, Bh, Bv]
# end

# function Cexp_old(x,y,z,v,args::WB_Model)
#   @unpack C, gNa, ENa, gK, EK, ϕ = args
#   n = ninf(v)
#   h = hinf(v)

#   dFdv3 = gNa*h*(-18*minf(v)*dminf(v)^2-9*minf(v)^2*d2minf(v)+6*dminf(v)^3*(ENa-v)+18*minf(v)*dminf(v)*d2minf(v)*(ENa-v)+3*minf(v)^2*d3minf(v)*(ENa-v))/C
#   dFdvdn2 = -12*gK*n^2/C
#   dFdn3 = 24*gK*n*(EK-v)/C
#   dFdv2dh = gNa*(-6*minf(v)^2*dminf(v)+6*minf(v)*dminf(v)^2*(ENa-v)+3*minf(v)^2*d2minf(v)*(ENa-v))/C

#   dNdv2dn = -1/d2τn(v,ϕ)
#   dNdv3 = d3ninf(v)/τn(v,ϕ)+3*d2ninf(v)/dτn(v,ϕ)+3*dninf(v)/d2τn(v,ϕ)

#   dHdv2dh = -1/d2τh(v,ϕ)
#   dHdv3 = d3hinf(v)/τh(v,ϕ)+3*d2hinf(v)/dτh(v,ϕ)+3*dhinf(v)/d2τh(v,ϕ)

#   Cv = dFdv3*x[3]*y[3]*z[3]+dFdvdn2*(x[3]*y[1]*z[1]+x[1]*y[3]*z[1]+x[1]*y[1]*z[3])+dFdn3*x[1]*y[1]*z[1]+dFdv2dh*(x[3]*y[3]*z[2]+x[3]*y[2]*z[3]+x[2]*y[3]*z[3])
#   Cn = dNdv2dn*(x[3]*y[3]*z[1]+x[3]*y[1]*z[3]+x[1]*y[3]*z[3])+dNdv3*x[3]*y[3]*z[3]
#   Ch = dHdv2dh*(x[3]*y[3]*z[2]+x[3]*y[2]*z[3]+x[2]*y[3]*z[3])+dHdv3*x[3]*y[3]*z[3]
#   return [Cn, Ch, Cv]
# end