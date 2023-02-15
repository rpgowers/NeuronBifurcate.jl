using Roots, ForwardDiff, Setfield, LinearAlgebra, NLsolve

vtrap(A, x, y) = ( abs(x) > eps() ? A*x/expm1(x/y) : A*(y-x/2+x^2/(12*y)-x^4/(720*y^3)) ) # mirrors the vtrap used in the Allen model

include("model_arguments.jl")
include("ml_functions.jl")
include("wb_functions.jl")

function soma_voltage(x,args)
  @unpack C, EL, gL = args
  v = x[end]
  return (gL*(EL-v)+Ia(x,args))/C
end

function vector_hessian(f, x)
  n = length(x)
  out = ForwardDiff.jacobian(x -> ForwardDiff.jacobian(f, x), x)
  return reshape(out, n, n, n)
end
export vector_hessian

function vector_hyper_hessian(f, x)
  n = length(x)
  H(X) =  ForwardDiff.jacobian(X -> ForwardDiff.jacobian(f, X), X)
  out = ForwardDiff.jacobian(x -> H(x), x)
  return reshape(out, n, n, n, n)
end
export vector_hyper_hessian

function I∞(v,args::Union{MLS_Param, WBS_Param})
  return soma_voltage(X∞(v,args),args)*args.C+args.Iext
end

function I∞(v, args::Union{MLDS_Param, WBDS_Param})
  @unpack C, ρ, xin, λ, gL, EL, Iext = args
  return soma_voltage(X∞(v,args),args)*C+ρ*gL*(EL-v)+Iext*exp(-xin/λ)
end

function bt(args::Union{MLS_Param, WBS_Param})
  @unpack C, dims = args
  F(x) = Ia(x, args)
  ∇F(x) = ForwardDiff.gradient(F, x)
  fbt(v) = 1+sum(τa(v,args).*∇F(X∞(v,args))[1:dims-1].*dA∞(v,args))/C  
  vbt = find_zeros(v-> fbt(v) , -80.0, 40.0)
  args_temp = @set args.gL = 0.0
  f2(v) = I∞(v,args_temp)
  gbt = ForwardDiff.derivative.(f2, vbt)
  Ibt = zeros(length(gbt))
  for i in eachindex(vbt)
    args_temp = setproperties(args, (gL=gbt[i], Iext=0.0))
    Ibt[i] = -I∞(vbt[i],args_temp)# /Iscale(args_temp)
  end
  return vbt, Ibt, gbt
end

function bt(args::Union{MLDS_Param, WBDS_Param})
  @unpack C, gL, τδ, dims = args
  F(x) = Ia(x, args)
  ∇F(x) = ForwardDiff.gradient(F, x)
  fbt(v) = C+0.5*τδ*(-gL+∇F(X∞(v,args))[end])+sum((0.5*τδ .+ τa(v,args)).*∇F(X∞(v,args))[1:dims-1].*dA∞(v,args))
  ρ(v) = (-gL+∇F(X∞(v,args))[end]+sum(∇F(X∞(v,args))[1:dims-1].*dA∞(v,args)))/gL
  vbt = find_zeros(v -> fbt(v), -80.0, 40.0)
  ρbt = ρ.(vbt)
  Ibt = zeros(length(vbt))
  for i in eachindex(Ibt)
    args_temp = setproperties(args, (ρ = ρbt[i], Iext=0.0))
    Ibt[i] = -I∞(vbt[i],args_temp)/Iscale(args_temp)
  end
  return vbt, Ibt, ρbt
end

export F
export soma_voltage
export I∞
export X∞
export bt
export hopf
export τa
export dA∞
export Ia

function Jacobian(x,args)
  f((x)) = F(x,args)
  J = ForwardDiff.jacobian(X -> f(X), x)
end

function vfps(args)
  vfp = find_zeros(v->I∞(v,args), -120.0, 50.0)
end
export vfps

function Iscale(args::Union{MLS_Param, WBS_Param})
  1
end

function Iscale(args::Union{MLDS_Param, WBDS_Param})
  @unpack xin, λ = args
  exp(-xin/λ)
end

function sn(args)
  args_temp = @set args.Iext = 0.0
  f(v) = I∞(v,args_temp)
  vsn = find_zeros(v->ForwardDiff.derivative(f,v), -100.0, 20.0)
  Isn = [-I∞(vsn[i],args_temp)/Iscale(args_temp) for i in eachindex(vsn)]
  return vsn, Isn
end
export sn

function cusp(args::Union{MLS_Param, WBS_Param})
  args_temp = setproperties(args, (gL=0.0, Iext=0.0))
  f(v) = I∞(v,args_temp)
  f1(v) = ForwardDiff.derivative(f,v)
  f2(v) = ForwardDiff.derivative(f1,v)

  vc = find_zeros(v->f2(v), -100.0, 20.0)
  gc = f1.(vc)
  Ic = zeros(length(gc))
  for i in eachindex(gc)
    args_temp = @set args.gL = gc[i]
    Ic[i] = -I∞.(vc[i],Ref(args_temp))/Iscale(args_temp)
  end
  return vc, Ic, gc
end

function cusp(args::Union{MLDS_Param, WBDS_Param})
  args_temp = setproperties(args, (ρ=0.0, Iext=0.0))
  f(v) = I∞(v,args_temp)
  f1(v) = ForwardDiff.derivative(f,v)
  f2(v) = ForwardDiff.derivative(f1,v)
  vc = find_zeros(v->f2(v), -100.0, 20.0)
  gc = f1.(vc)
  ρc = gc./args.gL
  Ic = zeros(length(ρc))
  for i in eachindex(ρc)
    args_temp = @set args.ρ = ρc[i]
    Ic[i] = -I∞.(vc[i],Ref(args_temp))/Iscale(args_temp)
  end
  return vc, Ic, ρc
end
export cusp

function btc(args::Union{MLS_Param, WBS_Param})
  @unpack dims = args
  vbtc, Ibtc, gbtc = cusp(args)
  F(x) = Ia(x, args)
  ∇F(x) = ForwardDiff.gradient(F, x)
  fbt(v) = sum(τa(v,args).*∇F(X∞(v,args))[1:dims-1].*dA∞(v,args))
  Cbtc = -fbt.(vbtc)
  return vbtc, Ibtc, gbtc, Cbtc
end

function btc(args::Union{MLDS_Param, WBDS_Param})
  @unpack C, gL, dims = args
  vbtc, Ibtc, ρbtc = cusp(args)
  F(x) = Ia(x, args)
  ∇F(x) = ForwardDiff.gradient(F, x)
  τbt(v) = -2*( C+sum( τa(v,args).*∇F(X∞(v,args))[1:dims-1].*dA∞(v,args) ) )/( -gL+∇F(X∞(v,args))[end]+ sum( ∇F(X∞(v,args))[1:dims-1].*dA∞(v,args) ) )
  τbtc = τbt.(vbtc)
  return vbtc, Ibtc, ρbtc, τbtc
end
export btc

z(ω, τδ) = sqrt(2+2*sqrt(1+ω^2*τδ^2))
# u(ω, τδ) = sign(ω)*sqrt(-2+2*sqrt(1+ω^2*τδ^2)) # we never actually want ω to be negative so do we need the sign?
u(ω, τδ) = sqrt(-2+2*sqrt(1+ω^2*τδ^2))


function hopf(args::Union{MLDS_Param, WBDS_Param}; v0 =-10.0, ω0=0.05, ωtol = 1e-3)
  @unpack dims, τδ, C, gL, ρ = args
  F(x) = Ia(x, args)
  ∇F(x) = ForwardDiff.gradient(F, x)

  function DS_hopftest!(F, (v, ω))
    F[1] = -gL+∇F(X∞(v,args))[end]-0.5*z(ω, τδ)*ρ*gL+sum( ∇F(X∞(v,args))[1:dims-1].*dA∞(v,args)./(1 .+ω^2 .*τa(v,args).^2 ) )
    F[2] = C+0.5*u(ω, τδ)*ρ*gL/ω+sum( ∇F(X∞(v,args))[1:dims-1].*dA∞(v,args).*τa(v,args)./(1 .+ω^2 .*τa(v,args).^2 ) )
    # F[1] = ∂f∂n(v,args)*da∞(v,An,Δn)*ψa(v,ϕ,An,Δn)^2/(ψa(v,ϕ,An,Δn)^2+ω^2)+∂f∂v(v,args)-0.5*z(ω, τδ)*ρ*gL/C
    # F[2] = -∂f∂n(v,args)*da∞(v,An,Δn)*ψa(v,ϕ,An,Δn)/(ψa(v,ϕ,An,Δn)^2+ω^2)-1-0.5*u(ω, τδ)*ρ*gL/(C*ω)
  end
  output = nlsolve(DS_hopftest!, [v0;ω0], autodiff = :forward, ftol=1e-8)
  vh = output.zero[1]
  ωh = output.zero[2]
  if ωh < ωtol
    vh = NaN
    Ih = NaN
  else
    args_temp = @set args.Iext = 0.0
    Ih = -I∞(vh,args_temp)
  end
  return vh, Ih, ωh
end
export hopf

function hopf_stability(v, ωh, args::Union{MLDS_Param, WBDS_Param})
  @unpack gL, ρ, C, τδ, dims = args
  χ = ρ*gL/C
  f(x) = Ia(x,args)/C
  ∇f(x) = ForwardDiff.gradient(f, x)

  #calculation of F
  qσ = 1.0
  γ = sqrt(1+1im*ωh*τδ)
  κ = 1+χ*τδ/(2*γ)+sum( ∇f(X∞(v,args))[1:dims-1].*dA∞(v,args).*τa(v,args)./(1 .+1im*ωh*τa(v,args)).^2 )
  q = vcat(dA∞(v,args)./(1 .+1im*ωh*τa(v,args)), [1.0])
  p = vcat(∇f(X∞(v,args))[1:dims-1].*τa(v,args)./(1 .+1im*ωh*τa(v,args)), [1.0])./κ
  Cvec = Cexp(q,q,conj.(q),v,args)
  F = conj(p)'Cvec
  # println(F)

  # calculation of G
  α = sum(∇f(X∞(v,args))[1:dims-1].*dA∞(v,args))+∇f(X∞(v,args))[end]-gL/C
  S = 1/(α-χ)
  J = S.*vcat(dA∞(v,args), [1.0])*vcat(∇f(X∞(v,args))[1:dims-1].*τa(v,args), [1.0])'+Diagonal( -vcat( τa(v,args), [0.0] ) )
  # display(J)
  Bu = Bexp(q, conj.(q), v, args)
  u = J*Bu
  B2 = Bexp(q, u, v, args)
  G = 2*conj(p)'*B2
  # println(G)
  
  # calculation of H
  Ω = 2*ωh
  Γ = sqrt(1+1im*Ω*τδ)
  α1 = sum( ∇f(X∞(v,args))[1:dims-1].*dA∞(v,args)./(1 .+ 1im*Ω*τa(v,args)) )+∇f(X∞(v,args))[end]-gL/C
  S1 = 1/(α1-1im*Ω-Γ*χ)
  By = Bexp(q,q,v,args)
  Ta = τa(v,args)./(1 .+ 1im*Ω*τa(v,args))
  K = S1.*vcat(Ta.*dA∞(v,args)./τa(v,args), [1.0])*conj(vcat(Ta.*∇f(X∞(v,args))[1:dims-1], [1.0]))'+Diagonal( -vcat( Ta, [0.0] ) )
  # display(K)
  y = -K*By
  B3 = Bexp(y, conj(q), v, args)
  H = conj(p)'*B3
  # println(H)
  return real(F-G+H)/(2*ωh) 
end
export hopf_stability