@testset "Map path" begin
@testset "Map path CLZ short" begin
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    path = [
        LLA(51.80665409406621, 10.335359255147063),
        LLA(51.806410, 10.335425),
        LLA(51.805925724425855, 10.335097770385982),
        LLA(51.80526071662234, 10.334263247682527),
        LLA(51.80494642844135, 10.334420156901205),
        LLA(51.8046968207022, 10.335259688289772),
        LLA(51.8043692602778, 10.335160446571944),
        LLA(51.80400106293405, 10.335266393838086)
    ]
    times = [now()+i*Second(30) for i in 1:length(path)]
    points = [EverySingleStreet.GPSPoint(p, ZonedDateTime(time, TimeZone("UTC"))) for (p, time) in zip(path, times)]
    gpxfile = EverySingleStreet.GPXFile("test", points)
    candidates = EverySingleStreet.map_path(city_map, gpxfile)
    @test length(candidates) == 1
    street_names = [c.way.name for c in candidates[1]]
    @test street_names == ["Windmühlenstraße","Adolph-Roemer-Straße","Adolph-Roemer-Straße","Adolph-Roemer-Straße","Graupenstraße","Graupenstraße","Graupenstraße","Schulstraße"]
end

@testset "Streetpath CLZ short" begin 
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    path = [
        LLA(51.80665409406621, 10.335359255147063),
        LLA(51.806410, 10.335425),
        LLA(51.805925724425855, 10.335097770385982),
        LLA(51.80526071662234, 10.334263247682527),
        LLA(51.80494642844135, 10.334420156901205),
        LLA(51.8046968207022, 10.335259688289772),
        LLA(51.8043692602778, 10.335160446571944),
        LLA(51.80400106293405, 10.335266393838086)
    ]
    times = [now()+i*Second(30) for i in 1:length(path)]
    points = [EverySingleStreet.GPSPoint(p, ZonedDateTime(time, TimeZone("UTC"))) for (p, time) in zip(path, times)]
    gpxfile = EverySingleStreet.GPXFile("test", points)

    candidates = EverySingleStreet.map_path(city_map, gpxfile)
    streetpaths = EverySingleStreet.calculate_streetpath("test", 1, candidates[1], city_map)
    streetpath = streetpaths[1]
    # check the segment properties
    for segment in streetpath.segments
        from = segment.from 
        to = segment.to
        @test from.way.name == to.way.name
        @test from.way_is_reverse == to.way_is_reverse
    end
end

@testset "Map matching walked parts" begin
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}())
    nt = EverySingleStreet.map_matching(joinpath(@__DIR__, "..", "data", "strava_hamburg.json"), city_map.ways, walked_parts, "tmp_local_map.json")
    @test nt.added_kms ≈ 0.0
    rm("tmp_local_map")
end

end