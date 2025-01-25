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
    return (g = g, nodes = map.nodes, kd_tree = kd_tree, osm_id_to_node_id = map.osm_id_to_node_id)
end

"""
    update_weights!(g::SimpleWeightedGraph, map::AbstractSimpleMap, walked_parts::WalkedParts; mul_non_walkable_road=2.5, mul_walked_road = 3.0)

Update the weights of the edges such that walked parts are "longer" in a shortest path calculation as well 
as ways that are not classified as `iswalkable_road`. 
"""
function update_weights!(g::SimpleWeightedGraph, map::AbstractSimpleMap, walked_parts::WalkedParts; mul_non_walkable_road=2.5, mul_walked_road = 3.0)
    for way in map.ways
        if !iswalkable_road(way)
            for i in 1:(length(way.nodes) - 1)
                u = map.osm_id_to_node_id[way.nodes[i].id]
                v = map.osm_id_to_node_id[way.nodes[i + 1].id]
                distance = euclidean_distance(way.nodes[i].lla, way.nodes[i+1].lla)
                distance *= mul_non_walkable_road
                add_edge!(g, u, v, distance)
            end
        end
        if haskey(walked_parts.ways, way.id)
            walked_way = walked_parts.ways[way.id]
            for i in 1:(length(walked_way.way.nodes) - 1)
                u = map.osm_id_to_node_id[way.nodes[i].id]
                v = map.osm_id_to_node_id[way.nodes[i + 1].id]
                c1 = node_on_way_to_candidate(i, walked_way.way)
                c2 = node_on_way_to_candidate(i+1, walked_way.way)
                intersections = Vector{Tuple{Float64, Float64}}()
                for part in walked_way.parts 
                    inters = intersection(part, (c1.λ,c2.λ))
                    isnothing(inters) && continue 
                    push!(intersections, inters)
                end
                d_way = euclidean_distance(walked_way.way.nodes[i].lla, walked_way.way.nodes[i+1].lla)
                d_inters = sum(p->p[2]-p[1], intersections; init=0)
                add_edge!(g, u, v,  d_way + d_inters * (mul_walked_road - 1))
            end
        end
    end
end

"""
    best_route(graph_data::Dict, src::LLA, dst::LLA)

Give a dict with keys "g", "nodes" and "kd_tree" compute the shortest path from `src` to `dst`.
The start and end point are first mapped to the graph and then a standard AStar is used.
Return a vector of `LLA` including the start and end point.
"""
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