@testset "GPX" begin 

@testset "Read gpx data" begin 
    path = EverySingleStreet.parse_gpx(joinpath(@__DIR__, "..", "data", "test.gpx"))
    @test path isa EverySingleStreet.GPXFile
    @test path.name == "test"
    @test path.gps_points[1].pos.lat ≈ 53.5578920 
    @test path.gps_points[1].pos.lon ≈ 9.9634340
end

end