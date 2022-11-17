using NeuronBifurcate, PyPlot

gin = 2.0:0.05:9.0
gσ = 2.0
ρ = gin./gσ .-1

τδ = [2.5, 5.0, 10.0, 20.0]

Isn = zeros(length(gin),2)
Ibt = zeros(length(τδ))
gbt = zeros(length(τδ))

vc, Ic, ρc = cusp(MLDS_Param())
gc = gσ.*(1 .+ρc)

println(Ic)
println(gc)

for j in eachindex(τδ)
  vbt, Ibt[j], ρbt = [bt(MLDS_Param(τδ=τδ[j]))[k][1] for k=1:3]
  gbt[j] = gσ*(1+ρbt)
  for i in eachindex(gin)
    args = MLDS_Param(gL=gσ, ρ=ρ[i], τδ=τδ[j])
    vout, Iout = sn(args)
    length(vout) > 1 ? Isn[i,:] .= Iout[1:2] : Isn[i,:] .= NaN
  end
end

println(Ibt)
println(gbt)

plot(Isn[:,1], gin, color="C0", label="SN")
plot(Isn[:,2], gin, color="C0")
plot(Ic[1], gc[1], "^", color="k", label="Cusp")

for j in eachindex(τδ)
  plot(Ibt[j], gbt[j], "o", color="C$(j)", alpha=0.5)
end

show()