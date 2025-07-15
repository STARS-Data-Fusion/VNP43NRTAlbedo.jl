using Test
using VNP43NRTAlbedo

@testset "VNP43NRTAlbedo basic tests" begin
    # Example test: check that the module loads and a key function exists
    @test isdefined(VNP43NRTAlbedo, :NRT_BRDF_all)
    # Add more tests here as needed
end
