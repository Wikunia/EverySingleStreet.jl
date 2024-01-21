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
    nt = EverySingleStreet.map_matching(joinpath(@__DIR__, "..", "data", "strava_hamburg.json"), city_map, walked_parts, "tmp_local_map.json")
    @test nt.added_kms ≈ 0.0
    rm("tmp_local_map.json")
end

@testset "Disconnected" begin 
    path = joinpath(@__DIR__, "..", "data", "Luebeck.json");
    EverySingleStreet.download("Lübeck, Germany", path);
    EverySingleStreet.filter_walkable_json!(path);
    _, city_map = EverySingleStreet.parse_no_graph_map(path, joinpath(@__DIR__, "..", "data", "luebeck_districts.geojson"));

    walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}())
    gps_points = EverySingleStreet.get_gps_points(joinpath(@__DIR__, "..", "data", "strava_luebeck.json"))
    gps_points = vcat(gps_points[1250:1280], gps_points[9680:9780])

    nt = EverySingleStreet.map_matching(gps_points, city_map, walked_parts, "tmp_local_map.json");
    rm("tmp_local_map.json")

    # Check that both parts exists even though there is no connection between them in the local map
    way_ids = nt.walked_parts.names["Große Altefähre"]
    ways = filter(p->first(p) in way_ids, walked_parts.ways)
    len = 0.0
    for (idx,way) in ways 
        for part in way.parts
            len += part[2]-part[1]
        end
    end
    @test len > 5

    way_ids = nt.walked_parts.names["Rehderbrücke"]
    ways = filter(p->first(p) in way_ids, walked_parts.ways)
    len = 0.0
    for (idx,way) in ways 
        for part in way.parts
            len += part[2]-part[1]
        end
    end
    @test len > 5

end

@testset "Altona walk" begin 
    _, altona_map = EverySingleStreet.parse_no_graph_map(joinpath(@__DIR__, "..", "data", "Altona.json"));
    altona_walk_path = joinpath(@__DIR__, "..", "data", "altona_walk.json");
    walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}());
    mm_data = EverySingleStreet.map_matching(altona_walk_path, altona_map, walked_parts, "tmp_local_map.json");
    rm("tmp_local_map.json")
    walked_parts = mm_data.walked_parts
    # Schillerstraße completely walked
    walked_way = walked_parts.ways[112450623]
    total_street_len = walked_way.way.meters
    walked_meters = sum(p[2]-p[1] for p in walked_way.parts)
    @test total_street_len ≈ walked_meters

    # Part in park is walked but should appear as it isn't a street
    @test !haskey(walked_parts.ways, 135511621)
end

end