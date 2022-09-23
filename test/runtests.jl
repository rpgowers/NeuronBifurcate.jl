using NeuronBifurcate
using Test

@testset "NeuronBifurcate.jl" begin
    # Write your tests here.
    args = MLS_Param()
    v = -50.0
    n = 0.1
    @test NeuronBifurcate.a∞(v, args.An, args.Δn) > 0
    @test NeuronBifurcate.ψa(v, args.ϕ, args.An, args.Δn) > 0
    @test NeuronBifurcate.ML_ncurrent((n,v),args) > -Inf
    @test NeuronBifurcate.MLS_voltage((n,v),args) > -Inf
    @test length(NeuronBifurcate.F((n,v), args)) == 2
    @test NeuronBifurcate.I∞(v,args) > -Inf
    @test length(NeuronBifurcate.vfps(args)) == 3
end
