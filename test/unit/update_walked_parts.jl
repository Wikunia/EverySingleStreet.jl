@testset "Add ends to street if nearly finished" begin 
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    path = [
        LLA(51.80513, 10.34662),
        LLA(51.80545, 10.34700),
        LLA(51.80570, 10.34769)
    ]
    times = [now()+i*Second(30) for i in 1:length(path)]
    points = [EverySingleStreet.GPSPoint(p, ZonedDateTime(time, TimeZone("UTC"))) for (p, time) in zip(path, times)]

    streetpaths = EverySingleStreet.map_matching(city_map, "test", points)
    wp = EverySingleStreet.streetpaths_to_walked_parts(streetpaths, city_map.ways)
    @test length(wp.ways) == 1
    @test haskey(wp.ways, 172200103)
    wway = wp.ways[172200103]
    # the street shouldn't be fully finished
    @test EverySingleStreet.total_length(wp) < wway.way.meters
    @test wway.parts[1][2] - wway.parts[1][1] ≈ EverySingleStreet.total_length(wp)

    prev_max = EverySingleStreet.get_preference("EXTEND_WALKED_WAY_UP_TO")
    EverySingleStreet.set_preferences!("EXTEND_WALKED_WAY_UP_TO" => 5)

    EverySingleStreet.extend_walked_parts!(wp)
    wway = wp.ways[172200103]
    # the street shouldn't be fully finished because of the preference of just adding up to 5m
    @test EverySingleStreet.total_length(wp) < wway.way.meters
    @test wway.parts[1][2] - wway.parts[1][1] ≈ EverySingleStreet.total_length(wp)

    EverySingleStreet.set_preferences!("EXTEND_WALKED_WAY_UP_TO" => 8)
    EverySingleStreet.extend_walked_parts!(wp)
    wway = wp.ways[172200103]
    # the way should be finished up til the end but not the start
    @test EverySingleStreet.total_length(wp) < wway.way.meters
    @test wway.parts[1][2] - wway.parts[1][1] ≈ EverySingleStreet.total_length(wp)
    @test wway.parts[1][2] ≈ wway.way.meters

    EverySingleStreet.set_preferences!("EXTEND_WALKED_WAY_UP_TO" => 12)
    EverySingleStreet.extend_walked_parts!(wp)
    wway = wp.ways[172200103]
    # the way should be completed
    @test EverySingleStreet.total_length(wp) ≈ wway.way.meters
    @test wway.parts[1][2] - wway.parts[1][1] ≈ EverySingleStreet.total_length(wp)
    @test isapprox(wway.parts[1][1], 0.0; atol=0.01)
    @test wway.parts[1][2] ≈ wway.way.meters

    EverySingleStreet.set_preferences!("EXTEND_WALKED_WAY_UP_TO" => prev_max)
end