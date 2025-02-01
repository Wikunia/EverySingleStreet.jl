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

    path = [
        LLA(51.811368,10.337099),
        LLA(51.810917, 10.336504),
        LLA(51.810548, 10.336252),
        LLA(51.809935, 10.336069),
    ]
    xdoc = EverySingleStreet.create_gpx_document(path)
    xroot = LightXML.root(xdoc)
    @test length(LightXML.get_elements_by_tagname(xroot, "trk")) == 1
    free(xdoc)


    times = [now()+i*Second(30) for i in 1:length(path)]
    points = [EverySingleStreet.GPSPoint(p, ZonedDateTime(time, TimeZone("UTC"))) for (p, time) in zip(path, times)]

    mapped = EverySingleStreet.map_matching(points, city_map, EverySingleStreet.WalkedParts(), "tmp_local_map.json");
    rm("tmp_local_map.json")

    best_points = EverySingleStreet.best_route(d, LLA(51.811368,10.337099), LLA(51.809935, 10.336069))
    before_update = EverySingleStreet.total_length(best_points)

    id1 = nt.osm_id_to_node_id[30337960]
    id2 = nt.osm_id_to_node_id[30337959]
    old_dist = nt.g.weights[id1, id2]
    EverySingleStreet.update_weights!(nt.g, city_map, mapped.walked_parts; mul_non_walkable_road=20, mul_walked_road=20)
   
    @test 20*old_dist ≈ nt.g.weights[id1, id2]

    d = Dict(
        "g" => nt.g,
        "nodes" => nt.nodes,
        "kd_tree" => nt.kd_tree,
    )
    best_points = EverySingleStreet.best_route(d, LLA(51.811368,10.337099), LLA(51.809935, 10.336069))
    @test before_update <= EverySingleStreet.total_length(best_points) <= 20*before_update
end