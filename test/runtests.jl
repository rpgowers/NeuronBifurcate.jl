using NeuronBifurcate
using Test

@testset "NeuronBifurcate.jl" begin
    # Write your tests here.
    args_mls = MLS_Param()
    v = -50.0
    n = 0.1
    @test a∞(v, args_mls.An, args_mls.Δn) > 0
    @test ψa(v, args_mls.ϕ, args_mls.An, args_mls.Δn) > 0
    @test ML_ncurrent((n,v),args_mls) > -Inf
    @test soma_voltage((n,v),args_mls) > -Inf
    @test length(F((n,v), args_mls)) == 2
    @test I∞(v,args_mls) > -Inf
    @test length(vfps(args_mls)) == 3
    @test length(sn(args_mls)[1]) == 2
    @test cusp(args_mls)[3][1] > 0.0
    @test bt(args_mls)[3][1] > 0.0
    @test hopf(args_mls)[3][2] > 0.0

    args_wbs = WBS_Param()
    v = -50.0
    n = 0.1
    h = 0.2

    @test soma_voltage((n,h,v),args_wbs) > -Inf
    @test length(F((n,h,v), args_wbs)) == 3
    @test I∞(v,args_wbs) > -Inf
    @test length(vfps(args_wbs)) == 3
    @test length(sn(args_wbs)[1]) == 2
    @test cusp(args_wbs)[3][1] > 0.0
    @test bt(args_wbs)[3][1] > 0.0

    args_mlds = MLDS_Param(ρ = 1.0)
    v = -50.0
    n = 0.1
    @test soma_voltage((n,h,v),args_mlds) > -Inf
    @test I∞(v,args_mlds) > -Inf
    @test length(vfps(args_mlds)) == 1
    @test length(sn(args_mlds)[1]) == 2
    @test cusp(args_mlds)[3][1] > 0.0
    @test bt(args_mlds)[3][1] > 0.0

end
