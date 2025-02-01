@testset "Candidates" begin

@testset "Single point several candidates" begin
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    p = EverySingleStreet.GPSPoint(LLA(51.80665409406621, 10.335359255147063), ZonedDateTime(now(), TimeZone("UTC")))
    candidates = EverySingleStreet.get_candidates(city_map, [p])[1]
    EverySingleStreet.points2geojson([candidates], "test.geojson")
    geojson = GeoJSON.read(read("test.geojson", String))
    @test length(geojson.geometry) == length(candidates)
    
    EverySingleStreet.create_connection_geojson([p], [candidates], "test.geojson")
    geojson = GeoJSON.read(read("test.geojson", String))
    @test length(geojson.geometry) == 1+2*length(candidates)
    rm("test.geojson")
    # check whether the following four streets have at least one candidate point
    test_street_names = ["Windmühlenstraße", "Kronenplatz", "Sorge", "Erzstraße"]
    for name in test_street_names
        found = false
        for candidate in candidates
            if candidate.way.name == name
                found = true
                break
            end
        end
        @test found
    end

    # check that lla of the candidate is less than 100m away from the initial point p
    max_dist = 0.0
    for candidate in candidates
        dist = euclidean_distance(p.pos, candidate.lla)
        max_dist = max(max_dist, dist)
        @test dist < 100
    end
    max_dist *= 0.75
    prev_max = EverySingleStreet.get_preference("CANDIDATES_MAXIMUM_DISTANCE")
    EverySingleStreet.set_preferences!("CANDIDATES_MAXIMUM_DISTANCE" => max_dist)
    candidates = EverySingleStreet.get_candidates(city_map, [p])[1]
    for candidate in candidates
        dist = euclidean_distance(p.pos, candidate.lla)
        @test dist < max_dist
    end
    # set back
    EverySingleStreet.set_preferences!("CANDIDATES_MAXIMUM_DISTANCE" => prev_max)
end

end