using Test

if !isdefined(@__MODULE__, :timelike_splittings_func)
    include("../core/timelike_splitting.jl")
end

@testset "nf-dependent timelike splitting kernels" begin
    x = 0.37

    @test timelike_splittings_func(x = x, type = "Pqg0", nf = 0) == 0.0
    @test timelike_splittings_func(x = x, type = "Pqg1", nf = 0) == 0.0
    @test timelike_splittings_func(x = x, type = "Pqg2", nf = 0) == 0.0

    _, pqq2_delta, pqq2_D0 = timelike_splittings_func(
        x = x,
        type = "Pqq2",
        nf = 5,
    )
    @test pqq2_delta ≈ (
        1295.625 - 5 * 173.935
        - 25 * (64 / 81) * (51 / 16 + 3 * z3 - 5 * z2)
    )
    @test pqq2_D0 ≈ 1174.898 - 5 * 183.187 - 25 * (64 / 81)

    _, pgg2_delta, pgg2_D0 = timelike_splittings_func(
        x = x,
        type = "Pgg2",
        nf = 5,
    )
    @test pgg2_delta ≈ 4425.451 - 5 * 528.719 + 25 * 6.4628
    @test pgg2_D0 ≈ 2643.521 - 5 * 412.172 - 25 * (16 / 9)

    coefficient = timelike_splitting_coefficient_func(
        y = x,
        loop_order = 2,
        nf = 5,
    )
    direct_pqg2 = timelike_splittings_func(x = x, type = "Pqg2", nf = 5)
    @test coefficient.Pqg ≈ direct_pqg2

    @test_throws DomainError timelike_splittings_func(
        x = x,
        type = "Pqg0",
        nf = -1,
    )
    @test_throws ArgumentError timelike_splittings_func(
        x = x,
        type = "not-a-kernel",
        nf = 5,
    )
end
