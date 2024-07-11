@testset "Drawing route" begin
    _, altona_map = EverySingleStreet.parse_no_graph_map(joinpath(@__DIR__, "..", "data", "Altona.json"));
    altona_walk_path = joinpath(@__DIR__, "..", "data", "altona_walk.json");
    walked_parts = EverySingleStreet.WalkedParts(Dict{String, Vector{Int}}(), Dict{Int, EverySingleStreet.WalkedWay}());
    mm_data = EverySingleStreet.map_matching(altona_walk_path, altona_map, walked_parts, "tmp_local_map.json");
    rm("tmp_local_map.json")
    walked_parts = mm_data.walked_parts

    gps_points = EverySingleStreet.get_gps_points(altona_walk_path)
    tmp_path = joinpath(@__DIR__, "..", "images", "tmp", "Altona_draw.png")
    EverySingleStreet.draw(walked_parts, gps_points, tmp_path)
    rm(tmp_path)
end