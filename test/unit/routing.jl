@testset "routing CLZ" begin 
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    nt = EverySingleStreet.convert_to_weighted_graph(city_map)
    @test nv(nt.g) > 20000
    id1 = nt.osm_id_to_node_id[2582004645]
    id2 = nt.osm_id_to_node_id[2582004646]
    @test nt.g.weights[id1, id2] > 10

    id1 = nt.osm_id_to_node_id[30223967]
    id2 = nt.osm_id_to_node_id[30223962]
    @test 30 <= nt.g.weights[id1, id2] <= 40

    id1 = nt.osm_id_to_node_id[30223967]
    id2 = nt.osm_id_to_node_id[30223962]
    old_dist = nt.g.weights[id1, id2]
    @test 30 <= nt.g.weights[id1, id2] <= 40

    EverySingleStreet.update_weights!(nt.g, city_map, EverySingleStreet.WalkedParts(); mul_non_walkable_road=5.0)
    @test nt.g.weights[id1, id2] ≈ 5*old_dist

    EverySingleStreet.update_weights!(nt.g, city_map, EverySingleStreet.WalkedParts(); mul_non_walkable_road=1.0)
    d = Dict(
        "g" => nt.g,
        "nodes" => nt.nodes,
        "kd_tree" => nt.kd_tree,
    )
    from = LLA(51.809676,10.335281)
    to = LLA(51.811626, 10.336118)
    best_points = EverySingleStreet.best_route(d, LLA(51.809676,10.335281), LLA(51.811626, 10.336118))
    @test best_points[1] == from
    @test best_points[end] == to
    @test 200 <= EverySingleStreet.total_length(best_points) <= 230

    EverySingleStreet.update_weights!(nt.g, city_map, EverySingleStreet.WalkedParts(); mul_non_walkable_road=5.0)
    d = Dict(
        "g" => nt.g,
        "nodes" => nt.nodes,
        "kd_tree" => nt.kd_tree,
    )
    best_points = EverySingleStreet.best_route(d, LLA(51.809676,10.335281), LLA(51.811626, 10.336118))
    @test best_points[1] == from
    @test best_points[end] == to
    @test 300 <= EverySingleStreet.total_length(best_points) <= 350
end