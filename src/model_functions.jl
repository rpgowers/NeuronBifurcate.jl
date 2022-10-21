using Roots, ForwardDiff, Setfield

include("model_arguments.jl")
include("ml_functions.jl")
include("wb_functions.jl")

export F
export soma_voltage
export I∞

function vfps(args)
  vfp = find_zeros(v->I∞(v,args), -120.0, 50.0)
end
export vfps

function Iscale(args::Union{MLS_Param, WBS_Param})
  1
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
  # idx = findfirst(gc .> 0) # choose the first positive conductance?
  Ic = zeros(length(gc))
  for i in eachindex(gc)
    args_temp = @set args.gL = gc[i]
    Ic[i] = -I∞.(vc[i],Ref(args_temp))/Iscale(args_temp)
  end
  return vc, Ic, gc
end
export cusp

export bt