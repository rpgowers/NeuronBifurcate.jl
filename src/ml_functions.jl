function a∞(v,A,Δ)
  return 1/( 1+exp( -(v-A)/Δ ) )
end
export a∞

function da∞(v,A,Δ)
  return exp( -(v-A)/Δ )/Δ/( 1+exp( -(v-A)/Δ ) )^2
end

function n∞(v,args::ML_Model)
  @unpack An, Δn = args
  return 1/( 1+exp( -(v-An)/Δn ) )
end
export n∞

function m∞(v,args::ML_Model)
  @unpack Am, Δm = args
  return 1/( 1+exp( -(v-Am)/Δm ) )
end
export m∞

function dn∞(v,args::ML_Model)
  @unpack An, Δn = args
  return exp( -(v-An)/Δn )/Δn/( 1+exp( -(v-An)/Δn ) )^2
end
export dn∞

function dm∞(v,args::ML_Model)
  @unpack Am, Δm = args
  return exp( -(v-Am)/Δm )/Δm/( 1+exp( -(v-Am)/Δm ) )^2
end
export dm∞

function ψa(v,args::ML_Model)
  @unpack An, Δn, ϕ = args
  return ϕ*cosh( (v-An)/(4*Δn) )
end
export ψa

function ML_ncurrent((n,v),args)
  # @unpack An, Δn = args
  return (n∞(v,args)-n)*ψa(v,args)
  # return (a∞(v,An,Δn)-n)*ψa(v,args)
end
export ML_ncurrent

function Ia((n,v),args)
  @unpack gs, Es, Am, Δm, gf, Ef = args
  # gs*n*(Es-v)+gf*a∞(v,Am,Δm)*(Ef-v)
  gs*n*(Es-v)+gf*m∞(v,args)*(Ef-v)
end

function X∞(v,args::ML_Model)
  # [a∞(v,args.An,args.Δn), v]
  [n∞(v,args), v]
end

function dA∞(v,args::ML_Model)
  # da∞(v,args.An,args.Δn)
  dn∞(v,args)
end

function τa(v,args::ML_Model)
  1/ψa(v,args)
end

function Fσ_dyn(x, args::ML_Model)
  [ML_ncurrent(x,args), soma_voltage(x,args)]
end

function ∂Ia∂n(v,args::ML_Model)
  @unpack gs, Es = args
  gs*(Es-v)
end

function ∂Ia∂v(v,args::ML_Model)
  @unpack gs, An, Δn, gf, Ef, Am, Δm = args
  # -gs*a∞(v,An,Δn)-gf*a∞(v,Am,Δm)+gf*da∞(v,Am,Δm)*(Ef-v)
  -gs*n∞(v,args)-gf*m∞(v,args)+gf*dm∞(v,args)*(Ef-v)
end

function hopf(args::MLS_Param)
  @unpack An, Δn, gL, C = args
  f(v) = -gL+∂Ia∂v(v,args)-C*ψa(v,args)
  vh = find_zeros(v->f(v), -80.0, 30.0)
  # ω2(v) = -ψa(v,args)*(ψa(v,args)+∂Ia∂n(v,args)*da∞(v,An,Δn)/C)
  ω2(v) = -ψa(v,args)*(ψa(v,args)+∂Ia∂n(v,args)*dn∞(v,args)/C)
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

function l1_coeff(p,q,B,C,ω,V,args)
  g20 = p'*B(q,q,V,args)
  g11 = p'*B(q,conj.(q),V,args)
  g21 = p'*C(q,q,conj.(q),V,args)
  real(1im*g20*g11+ω*g21)/(2*ω^2)
end

function hopf_stability(vh, ωh, args::MLS_Param)
  J = Jacobian(X∞(vh,args), args)
  q = eigvecs(J)[:,2]
  p = eigvecs(Array(J'))[:,1]
  k = p'*q
  pn = -p./k 
  l1_coeff(pn,q,Bexp,Cexp,ωh,vh,args)
end