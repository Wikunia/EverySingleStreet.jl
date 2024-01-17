@testset "get geojson file" begin
    districts = EverySingleStreet.get_districts(joinpath(@__DIR__, "..", "data", "hamburg_districts.geojson")); 
    @test length(districts) == 104
    district_names = [String(district.name) for district in districts]
    @test "Eimsbüttel" in district_names
    @test "Winterhude" in district_names
    @test "Bergedorf" in district_names
end

@testset "get district from node" begin
    nodes = [
        EverySingleStreet.Node(1, 53.57355, 9.95301) # in Eimsbüttel
        EverySingleStreet.Node(2, 53.5684, 9.9864) # in Rotherbaum
        EverySingleStreet.Node(3, 53.56155, 10.02219) # in Hohenfelde
    ]
    district_names = EverySingleStreet.map_nodes_to_district(nodes, joinpath(@__DIR__, "..", "data", "hamburg_districts.geojson"))
    @test district_names == [:Eimsbüttel, :Rotherbaum, :Hohenfelde]
end

@testset "Lübeck walked districts" begin 
    path = joinpath(@__DIR__, "..", "data", "Luebeck.json");
    EverySingleStreet.download("Lübeck, Germany", path);
    EverySingleStreet.filter_walkable_json!(path);
    _, city_map = EverySingleStreet.parse_no_graph_map(path, joinpath(@__DIR__, "..", "data", "luebeck_districts.geojson"));

    walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}())
    nt = EverySingleStreet.map_matching(joinpath(@__DIR__, "..", "data", "strava_luebeck.json"), city_map, walked_parts, "tmp_local_map.json");
    rm("tmp_local_map.json")
    @test nt.added_kms > 5
    district_percentages = EverySingleStreet.get_walked_district_perc(city_map, collect(values(nt.walked_parts.ways)))
    @test haskey(district_percentages, :Innenstadt)
    @test haskey(district_percentages, Symbol("Sankt Lorenz Süd"))

    # centroid of filtered walked parts 
    walked_parts = nt.walked_parts
    way_ids = nt.walked_parts.names["Große Altefähre"]
    ways = filter(p->first(p) in way_ids, walked_parts.ways)
    walked_parts = EverySingleStreet.WalkedParts(walked_parts.names, ways)
    centroid = EverySingleStreet.get_centroid(walked_parts)
    dist = euclidean_distance(centroid, LLA(53.873281, 10.688440))
    @test dist < 100
end