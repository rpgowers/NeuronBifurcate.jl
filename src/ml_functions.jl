function a∞(v,A,Δ)
  return 1/( 1+exp( -(v-A)/Δ ) )
end
export a∞

function da∞(v,A,Δ)
  return exp( -(v-A)/Δ )/Δ/( 1+exp( -(v-A)/Δ ) )^2
end

# function d2a∞(v,A,Δ) # this is just to compare against previous functions, delete when no longer needed
#   return 2*exp( -(v-A)/Δ )^2/Δ^2/( 1+exp( -(v-A)/Δ ) )^3-exp( -(v-A)/Δ )/Δ^2/( 1+exp( -(v-A)/Δ ) )^2
# end

# function d3a∞(v,A,Δ) # this is just to compare against previous functions, delete when no longer needed
#   val = 6*exp( -3*(v-A)/Δ )/( 1+exp( -(v-A)/Δ ) )^4-6*exp( -2*(v-A)/Δ )/( 1+exp( -(v-A)/Δ ) )^3+exp( -(v-A)/Δ )/( 1+exp( -(v-A)/Δ ) )^2
#   return val/Δ^3
# end

function ψa(v,ϕ,A,Δ)
  return ϕ*cosh( (v-A)/(4*Δ) )
end
export ψa

# function dψa(v,ϕ,A,Δ) # this is just to compare against previous functions, delete when no longer needed
#   return ϕ*sinh( (v-A)/(4*Δ) )/(4*Δ)
# end

# function d2ψa(v,ϕ,A,Δ) # this is just to compare against previous functions, delete when no longer needed
#   return ϕ*cosh( (v-A)/(4*Δ) )/(16*Δ^2)
# end

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

function F(x, args::ML_Model) # note that this does not currently include Iext in soma_voltage
  [ML_ncurrent(x,args), soma_voltage(x,args)]
end

function ∂Ia∂n(v,args::ML_Model)
  @unpack gs, Es = args
  gs*(Es-v)
end

function ∂Ia∂v(v,args::ML_Model)
  @unpack gs, An, Δn, gf, Ef, Am, Δm = args
  -gs*a∞(v,An,Δn)-gf*a∞(v,Am,Δm)+gf*da∞(v,Am,Δm)*(Ef-v)
end

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

# function Bexp_old(x,y,v,args::ML_Model) # this is the old function for finding B, delete when no longer needed
#   @unpack C, gs, gf, Ef, Am, Δm, An, Δn, ϕ = args
#   Bn = (2*dψa(v,ϕ,An,Δn)*da∞(v,An,Δn)+d2a∞(v,An,Δn)*ψa(v,ϕ,An,Δn))*x[2]*y[2]-dψa(v,ϕ,An,Δn)*(x[1]*y[2]+x[2]*y[1])
#   Bv = gf*(d2a∞(v,Am,Δm)*(Ef-v)-2*da∞(v,Am,Δm))/C*x[2]*y[2]-gs/C*(x[1]*y[2]+x[2]*y[1])
#   return [Bn, Bv]
# end
# export Bexp_old

function Bexp(x,y,v,args)
  @unpack dims = args
  f(X) = F(X,args)
  H = vector_hessian(f, X∞(v,args))
  return [sum(H[k,i,j]*x[i]*y[j] for i in eachindex(x), j in eachindex(y)) for k=1:dims]
end
export Bexp

# function Cexp_old(x,y,z,v,args::ML_Model)
#   @unpack C, gs, gf, Ef, Am, Δm, An, Δn, ϕ = args
#   Cn1 = (+3*d2ψa(v,ϕ,An,Δn)*da∞(v,An,Δn)+3*dψa(v,ϕ,An,Δn)*d2a∞(v,An,Δn)+d3a∞(v,An,Δn)*ψa(v,ϕ,An,Δn))*x[2]*y[2]*z[2]
#   Cn2 = -d2ψa(v,ϕ,An,Δn)*(x[2]*y[2]*z[1]+x[2]*y[1]*z[2]+x[1]*y[2]*z[2])
#   # the older values below seem to have been slightly incorrect?
#   # Cn1 = (-3*d2ψa(v,ϕ,An,Δn)*da∞(v,An,Δn)-3*dψa(v,ϕ,An,Δn)*d2a∞(v,An,Δn)+d3a∞(v,An,Δn)*ψa(v,ϕ,An,Δn))*x[2]*y[2]*z[2]
#   # Cn2 = -3*d2ψa(v,ϕ,An,Δn)*(x[2]*y[2]*z[1]+x[2]*y[1]*z[2]+x[1]*y[2]*z[2])
#   Cv = gf*(d3a∞(v,Am,Δm)*(Ef-v)-3*d2a∞(v,Am,Δm))/C*x[2]*y[2]*z[2]
#   return [Cn1+Cn2, Cv]
# end
# export Cexp_old

function Cexp(x,y,z,v,args)
  @unpack dims = args
  f(X) = F(X,args)
  V = vector_hyper_hessian(f, X∞(v,args))
  return [sum(V[i,j,k,l]*x[j]*y[k]*z[l] for j in eachindex(x), k in eachindex(y), l in eachindex(z) ) for i=1:dims ]
end
export Cexp

function l1_coeff(p,q,B,C,ω,V,args)
  g20 = p'*B(q,q,V,args)
  g11 = p'*B(q,conj.(q),V,args)
  g21 = p'*C(q,q,conj.(q),V,args)
  real(1im*g20*g11+ω*g21)/(2*ω^2)
end

function hopf_stability((nh, vh), ωh, args::MLS_Param)
  J = Jacobian([nh,vh],args)
  q = eigvecs(J)[:,2]
  p = eigvecs(Array(J'))[:,1]
  k = p'*q
  pn = -p./k 
  l1_coeff(pn,q,Bexp,Cexp,ωh,vh,args)
end

function ∂f∂v(v,args::Union{ML_Model})
  @unpack C, gL, gs, An, Δn, gf, Am, Δm, Ef = args
  (-gL-gs*a∞(v,An,Δn)-gf*a∞(v,Am,Δm)+gf*da∞(v,Am,Δm)*(Ef-v) )/C
end

function ∂f∂n(v,args::Union{ML_Model})
  @unpack gs, Es, C = args
  gs*(Es-v)/C
end

function hopf_stability((nfp,vfp),ωh,args::MLDS_Param)
  @unpack gL, ρ, C, τδ, ϕ, An, Δn = args
  τa = 1/ψa(vfp,ϕ,An,Δn)
  qσ = 1
  qa = qσ*da∞(vfp,An,Δn)/(1+1im*ωh*τa)
  κ = τa*∂f∂n(vfp,args)*da∞(vfp,An,Δn)/(1+1im*ωh*τa)^2+1+τδ*ρ*gL/(2*C*sqrt(1+1im*ωh*τδ))
  pσ = 1/κ
  pa = pσ*τa*∂f∂n(vfp,args)/(1+1im*ωh*τa)
  Ca, Cσ = Cexp([qa,qσ],[qa,qσ],[conj(qa),qσ],vfp,args)
  F = pa*Ca+pσ*Cσ

  Bau, Bσu = Bexp([qa,qσ],[conj(qa),qσ],vfp,args)
  a1 = -ψa(vfp,ϕ,An,Δn)
  b1 = da∞(vfp,An,Δn)*ψa(vfp,ϕ,An,Δn)
  c1 = ∂f∂n(vfp, args)
  a2 = ∂f∂v(vfp, args)
  x = ρ*gL/C
  J11 = ( a2-x )/(-b1*c1+a1*( a2-x ) )
  ϕ3θM2 = 1/(-b1*c1+a1*(a2-x) )
  J12 = -b1*ϕ3θM2
  J21 = -c1*ϕ3θM2
  J22 = a1*ϕ3θM2
  ua = J11*Bau+J12*Bσu
  uσ = J21*Bau+J22*Bσu
  Ba2, Bσ2 = Bexp([qa,qσ],[ua,uσ],vfp,args)
  G = 2*(pa*Ba2+pσ*Bσ2)

  Bay, Bσy = Bexp([qa,qσ],[qa,qσ],vfp,args)
  ω = 2*ωh
  γ = sqrt(1+1im*ω*τδ)
  Den = -b1*c1*γ+(a1-1im*ω)*( (a2-1im*ω)*γ-x*(1+1im*ω*τδ) )
  K11 = ( (a2-1im*ω)*γ-x*(1+1im*ω*τδ) )/Den
  ϕ3θM2 = γ/Den
  K12 = -b1*ϕ3θM2
  K21 = -c1*ϕ3θM2
  K22 = (a1-1im*ω)*ϕ3θM2
  ya = -(K11*Bay+K12*Bσy)
  yσ = -(K21*Bay+K22*Bσy)
  Ba3, Bσ3 = Bexp([conj(qa), qσ],[ya,yσ],vfp,args)
  H = pa*Ba3+pσ*Bσ3
  real(F-G+H)/(2*ωh)
end

export hopf_stability