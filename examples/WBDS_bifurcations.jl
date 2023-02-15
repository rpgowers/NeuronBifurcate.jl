using NeuronBifurcate, PyPlot

function hopf_finder(v0, ω0, args; ωtol=1e-3, vtol=1e-2)
  hstop = 0
  vhl, Ihl, ωhl = hopf(args; v0=v0[1], ω0=ω0[1], ωtol = ωtol)
  vhu, Ihu, ωhu = hopf(args; v0=v0[2], ω0=ω0[2], ωtol = ωtol)
  if abs(vhu-vhl) < vtol
    Ihl = NaN
    Ihu = NaN
    hstop = 1
  end
  return [Ihl, Ihu], [ωhl, ωhu], [vhl, vhu], hstop
end

function hopfstab_finder(vh, ωh, Ih,args)
  Sh = [0.0, 0.0]
  if Ih[1] > -Inf
    Sh[1] = hopf_stability(vh[1],ωh[1],args)
  end
  if Ih[2] > -Inf
    Sh[2] = hopf_stability(vh[2],ωh[2],args)
  end
  return Sh
end

gin = 0.1:0.01:1.5
gσ = 0.1
ρ = gin./gσ .-1

τδ = [1.0, 2.5, 5.0, 10.0]

Isn = zeros(length(gin),2)
Ih = zeros(length(τδ), length(gin), 2)
Sh = zeros(length(τδ), length(gin), 2)
ωh = zeros(length(τδ), length(gin), 2)
vbt = zeros(length(τδ))
Ibt = zeros(length(τδ))
gbt = zeros(length(τδ))

vc, Ic, ρc = cusp(WBDS_Param())
gc = gσ.*(1 .+ρc)


hstop = 0
for j in eachindex(τδ)
  global hstop = 0
  vbt[j], Ibt[j], ρbt = [bt(WBDS_Param(τδ=τδ[j]))[k][1] for k=1:3]
  gbt[j] = gσ*(1+ρbt)
  for i in eachindex(gin)
    args = WBDS_Param(gL=gσ, ρ=ρ[i], τδ=τδ[j])
    vout, Iout = sn(args)
    length(vout) > 1 ? Isn[i,:] .= Iout[1:2] : Isn[i,:] .= NaN
    if hstop == 0
      # if gin[i] > gbt[j]
        Ih[j,i,:], ωh[j,i,:], vh, hstop = hopf_finder([vbt[j], -30.0], [1.5e-2, 1.7], args; ωtol=1e-2, vtol=1e-2)
        Sh[j,i,:] .= hopfstab_finder([vh[1], vh[2]], ωh[j,i,:], Ih[j,i,:],args)
      # else
      #   Ih[j,i,1] = NaN
      #   vh, Ih[j,i,2], ωh[j,i,2] = hopf(args; v0=-30.0, ω0=1.7, ωtol = 0.0)
        # Sh[j,i,2] = 1.0
      # end
    else
      Ih[j,i,:] .= NaN
    end
  end
end


plot(Isn[:,1], gin, color="C0", label="SN")
plot(Isn[:,2], gin, color="C0")
plot(Ic[1], gc[1], "^", color="k", label="Cusp")

for j in eachindex(τδ)
  if j==1
	  plot(Ih[j,:,1][findall(x->x>0, Sh[j,:,1])], gin[findall(x->x>0, Sh[j,:,1])], ":", color="C$(j)", alpha=0.5, label="Sub-Hopf")
	  plot(Ih[j,:,1][findall(x->x<0, Sh[j,:,1])], gin[findall(x->x<0, Sh[j,:,1])], "--", color="C$(j)", alpha=0.5, label="Sup-Hopf")
	else
	  plot(Ih[j,:,1][findall(x->x>0, Sh[j,:,1])], gin[findall(x->x>0, Sh[j,:,1])], ":", color="C$(j)", alpha=0.5)
	  plot(Ih[j,:,1][findall(x->x<0, Sh[j,:,1])], gin[findall(x->x<0, Sh[j,:,1])], "--", color="C$(j)", alpha=0.5)
	end

	plot(Ih[j,:,2][findall(x->x>0, Sh[j,:,2])], gin[findall(x->x>0, Sh[j,:,2])], ":", alpha=0.5, color="C$(j)")
	plot(Ih[j,:,2][findall(x->x<0, Sh[j,:,2])], gin[findall(x->x<0, Sh[j,:,2])], "--", alpha=0.5, color="C$(j)")
	plot(Ibt[j], gbt[j], "o", color="C$(j)", alpha=0.5, label="BT, \$\\tau_\\delta\$ = $(τδ[j]) ms")
end
legend(frameon=false)
show()