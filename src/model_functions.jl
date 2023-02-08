using Roots, ForwardDiff, Setfield, LinearAlgebra

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
  # out = ForwardDiff.jacobian(x -> Matrix(transpose(ForwardDiff.jacobian(f, x))), x)
  # return reshape(Matrix(transpose(out)), n, n, n)
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

function I∞(v, args::Union{MLDS_Param})
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

function bt(args::Union{MLDS_Param})
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
export bt
export hopf

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

function Iscale(args::Union{MLDS_Param})
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

function cusp(args::Union{MLDS_Param})
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
export btc