function fax((vσ, vδ), args::Union{MLMDS_Param, MLADS_Param}) # discretised axial current derivative
  @unpack ρ, gL, λ, C, L, M = args
  dx = L/M
  ρ*gL*λ*(vδ-vσ)/(dx*C)
end

function fδ((vm,v0,vp),args::Union{MLMDS_Param, MLADS_Param})
  @unpack M, L, λ, τδ, EL = args
  dx = L/M
  Λ = λ^2/dx^2
  (EL-v0+Λ*(vp-2*v0+vm))/τδ
end

function fseal((vm, v0), args::Union{MLMDS_Param, MLADS_Param})
  @unpack M, L, λ, τδ, EL = args
  dx = L/M
  Λ = λ^2/dx^2
  (EL-v0+2*Λ*(vm-v0))/τδ
end

function discrete_rho(gin, args::Union{MLMDS_Param})
  @unpack λ, L, M, gL = args
  Δx = L/M
  Λ = λ^2/Δx^2
  μp = ((1+2*Λ)+sqrt(1+4*Λ))/(2*Λ)
  μm = ((1+2*Λ)-sqrt(1+4*Λ))/(2*Λ)
  D = (μp^(M+1)-μp^(M-1))-(μm^(M+1)-μm^(M-1))
  F = ( (1-μm)*(μp^(M+1)-μp^(M-1))-(1-μp)*(μm^(M+1)-μm^(M-1)) )/D
  ρ = (gin-gL)/(gL*F*sqrt(Λ))
end
export discrete_rho

function discrete_gin(args::Union{MLMDS_Param})
  @unpack λ, L, M, gL, ρ = args
  Δx = L/M
  Λ = λ^2/Δx^2
  μp = ((1+2*Λ)+sqrt(1+4*Λ))/(2*Λ)
  μm = ((1+2*Λ)-sqrt(1+4*Λ))/(2*Λ)
  D = (μp^(M+1)-μp^(M-1))-(μm^(M+1)-μm^(M-1))
  F = ( (1-μm)*(μp^(M+1)-μp^(M-1))-(1-μp)*(μm^(M+1)-μm^(M-1)) )/D
  gin = gL*(1+ρ*F*sqrt(Λ))
end
export discrete_gin