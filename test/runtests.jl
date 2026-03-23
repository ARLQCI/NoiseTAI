using Test
using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports
using NoiseTAI

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(NoiseTAI) === nothing
    @test check_no_stale_explicit_imports(NoiseTAI) === nothing
end

include("test_noise_TAI.jl")
