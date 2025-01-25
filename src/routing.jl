"""
    convert_to_weighted_graph(map::NoGraphMap)

Create a [`SimpleWeightedGraph`](@ref) from an [`AbstractSimpleMap`](@ref). Use the euclidean_distance between two points on a way as the weight.
"""
function convert_to_weighted_graph(map::AbstractSimpleMap)
    g = SimpleWeightedGraph(length(map.nodes))
    for way in map.ways
        for i in 1:(length(way.nodes) - 1)
            u = map.osm_id_to_node_id[way.nodes[i].id]
            v = map.osm_id_to_node_id[way.nodes[i + 1].id]
            distance = euclidean_distance(way.nodes[i].lla, way.nodes[i+1].lla)
            add_edge!(g, u, v, distance)
        end
    end
    origin_lla = get_centroid(map.nodes)
    kd_tree = get_kd_tree_from_points(map.nodes, origin_lla)
    return (g = g, nodes = map.nodes, kd_tree = kd_tree)
end

function best_route(graph_data::Dict, src::LLA, dst::LLA)
    g = graph_data["g"]
    nodes = graph_data["nodes"]
    kdtree = graph_data["kd_tree"]
    origin_lla = get_centroid(nodes)
    trans = ENUfromLLA(origin_lla, wgs84) 
    src_id, _ = nn(kdtree, Point2{Float64}(getxy(src, trans)))
    dst_id, _ = nn(kdtree, Point2{Float64}(getxy(dst, trans)))
    astar_result = Graphs.Experimental.ShortestPaths.shortest_paths(g, src_id, dst_id, Graphs.Experimental.ShortestPaths.AStar())

    lla_points = [src]
    for id in astar_result.path
        push!(lla_points, nodes[id].lla)
    end
    push!(lla_points, dst)
    return lla_points
end