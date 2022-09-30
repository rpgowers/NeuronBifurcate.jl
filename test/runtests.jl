using NeuronBifurcate
using Test

@testset "NeuronBifurcate.jl" begin
    # Write your tests here.
    args = MLS_Param()
    v = -50.0
    n = 0.1
    @test a∞(v, args.An, args.Δn) > 0
    @test ψa(v, args.ϕ, args.An, args.Δn) > 0
    @test ML_ncurrent((n,v),args) > -Inf
    @test MLS_voltage((n,v),args) > -Inf
    @test length(F((n,v), args)) == 2
    @test I∞(v,args) > -Inf
    @test length(vfps(args)) == 3
    @test length(sn(args)[1]) == 2
end
