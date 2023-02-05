@testset "Candidates" begin

@testset "Single point several candidates" begin
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    p = EverySingleStreet.GPSPoint(LLA(51.80665409406621, 10.335359255147063), ZonedDateTime(now(), TimeZone("UTC")))
    candidates = EverySingleStreet.get_candidates(city_map, [p])[1]
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
    for candidate in candidates
        @test euclidean_distance(p.pos, candidate.lla) < 100
    end
end

end