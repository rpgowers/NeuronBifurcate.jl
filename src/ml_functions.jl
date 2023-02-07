function a∞(v,A,Δ)
  return 1/( 1+exp( -(v-A)/Δ ) )
end
export a∞

function da∞(v,A,Δ)
  return exp( -(v-A)/Δ )/Δ/( 1+exp( -(v-A)/Δ ) )^2
end

function d2a∞(v,A,Δ) # this is just to compare against previous functions, delete when no longer needed
  return 2*exp( -(v-A)/Δ )^2/Δ^2/( 1+exp( -(v-A)/Δ ) )^3-exp( -(v-A)/Δ )/Δ^2/( 1+exp( -(v-A)/Δ ) )^2
end

function ψa(v,ϕ,A,Δ)
  return ϕ*cosh( (v-A)/(4*Δ) )
end
export ψa

function dψa(v,ϕ,A,Δ) # this is just to compare against previous functions, delete when no longer needed
  return ϕ*sinh( (v-A)/(4*Δ) )/(4*Δ)
end

function d2ψa(v,ϕ,A,Δ) # this is just to compare against previous functions, delete when no longer needed
  return ϕ*cosh( (v-A)/(4*Δ) )/(16*Δ^2)
end

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

function Bexp_old(x,y,v,args::ML_Model) # this is the old function for finding B, delete when no longer needed
  @unpack C, gs, gf, Ef, Am, Δm, An, Δn, ϕ = args
  Bn = (2*dψa(v,ϕ,An,Δn)*da∞(v,An,Δn)+d2a∞(v,An,Δn)*ψa(v,ϕ,An,Δn))*x[2]*y[2]-dψa(v,ϕ,An,Δn)*(x[1]*y[2]+x[2]*y[1])
  Bv = gf*(d2a∞(v,Am,Δm)*(Ef-v)-2*da∞(v,Am,Δm))/C*x[2]*y[2]-gs/C*(x[1]*y[2]+x[2]*y[1])
  return [Bn, Bv]
end
export Bexp_old

function Bexp(x,y,v,args::ML_Model)
  @unpack dims = args
  f(X) = F(X,args)
  H = vector_hessian(f, X∞(v,args))
  return [sum(H[i,j,k]*x[i]*y[j] for i in eachindex(x), j in eachindex(y)) for k=1:dims]
end
export Bexp

function Cexp(x,y,z,v,args::ML_Model)
  @unpack C, gs, gf, Ef, Am, Δm, An, Δn, ϕ = args
  Cn1 = (-3*d2ψa(v,ϕ,An,Δn)*da∞(v,An,Δn)-3*dψa(v,ϕ,An,Δn)*d2a∞(v,An,Δn)+d3a∞(v,An,Δn)*ψa(v,ϕ,An,Δn))*x[2]*y[2]*z[2]
  Cn2 = -3*d2ψa(v,ϕ,An,Δn)*(x[2]*y[2]*z[1]+x[2]*y[1]*z[2]+x[1]*y[2]*z[2])
  Cv = gf*(d3a∞(v,Am,Δm)*(Ef-v)-3*d2a∞(v,Am,Δm))/C*x[2]*y[2]*z[2]
  return [Cn1+Cn2, Cv]
end

function l1_coeff(p,q,B,C,ω,V,args)
  g20 = p'*B(q,q,V,args)
  g11 = p'*B(q,conj.(q),V,args)
  g21 = p'*C(q,q,conj.(q),V,args)
  real(1im*g20*g11+ω*g21)/(2*ω^2)
end

function hopf_stability(args::MLS_Param)
  J = Jacobian([nh,vh],args)
  q = eigvecs(J)[:,2]
  p = eigvecs(Array(J'))[:,1]
  k = p'*q
  pn = -p./k 
  l1_coeff(pn,q,Bexp,Cexp,ωh,vh,args)
end