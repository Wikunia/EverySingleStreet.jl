"""
    bounded_dijkstra(graph, dist_mat, source, distance)

Compute the shortest distance and the parents on how to get there from the specified source up to distance away.
Return distances as a dictionary pointing from destination to shortest distance as well as parents point from destination on how to get there.
"""
function bounded_dijkstra(graph, dist_mat, source, distance)
    distances = Dict{Int32, Float64}()
    parents = Dict{Int32,Int32}()
    distances[source] = 0

    seen = Set{Int32}()
    queue = PriorityQueue{Int32, Float64}()
    enqueue!(queue, source => distances[source])

    while !isempty(queue)
        node = dequeue!(queue)
        push!(seen, node)

        for neighbor in neighbors(graph, node)
            if get(distances, neighbor, Inf) > get(distances, node, Inf) + dist_mat[node, neighbor] && get(distances, node, Inf) <= distance
                distances[neighbor] = get(distances, node, Inf) + dist_mat[node, neighbor]
                parents[neighbor] = node
                if !(neighbor in seen)
                    if haskey(queue, neighbor)
                        queue[neighbor] = distances[neighbor]
                    else
                        enqueue!(queue, neighbor => distances[neighbor])
                    end
                end
            end
        end
    end

    return distances, parents
end

"""
    get_dist_mat_walkable(map, osm_graph, dist_mat::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti}

Set the distance for non walkable road ids to the highest value such that they will not be used. 
"""
function get_dist_mat_walkable(map, osm_graph, dist_mat::SparseMatrixCSC{Tv, Ti}) where {Tv, Ti}
    dist_mat_walkable = deepcopy(dist_mat)
    for way in map.ways 
        iswalkable_road(way) && continue
        nodes = way.nodes
        for (i, j) in zip(1:length(nodes)-1,2:length(nodes)) 
            idx1 = get(osm_graph.node_to_index, nodes[i].id, nothing)
            idx2 = get(osm_graph.node_to_index, nodes[j].id, nothing)
            if isnothing(idx1) || isnothing(idx2)
                continue
            end 
            dist_mat_walkable[idx1, idx2] = typemax(Tv)
            dist_mat_walkable[idx2, idx1] = typemax(Tv)
        end
    end
    return dist_mat_walkable
end

"""
    bounded_all_shortest_paths(map::AbstractSimpleMap, osm_graph, distance)

Call [`bounded_dijkstra`](@ref) for all nodes and return a `BoundedAllShortestPaths` object.
"""
function bounded_all_shortest_paths(map::AbstractSimpleMap, osm_graph, distance)
    g = osm_graph.graph
    dist_mat = osm_graph.weights
    nodeid_to_local = map.osm_id_to_node_id
    walkable_road_nodes = map.walkable_road_nodes
    dist_mat_walkable = get_dist_mat_walkable(map, osm_graph, dist_mat)
   
    all_parents = Vector{Dict{Int32, Int32}}(undef, nv(g))
    all_parents_walable = Vector{Dict{Int32, Int32}}(undef, nv(g))
    @showprogress for i in 1:nv(g)
        if !walkable_road_nodes[nodeid_to_local[osm_graph.index_to_node[i]]]
            # all_distances[i] = Dict{Int32, Float64}()
            # all_parents[i] = Dict{Int32, Int32}()
            continue
        end
        distances, parents = bounded_dijkstra(g, dist_mat, i, distance)
        distances, parents_walkable = bounded_dijkstra(g, dist_mat_walkable, i, distance)
        all_parents[i] = parents
        all_parents_walable[i] = parents_walkable
    end
    return BoundedAllShortestPaths(g, dist_mat, dist_mat_walkable, all_parents, all_parents_walable, distance)
end

"""
    get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to; only_walkable_road=false)

Return the same output as `a_star(g, from, to, dist_mat)` but use the cache from `BoundedAllShortestPaths`
"""
function get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to; only_walkable_road=false)
    g = bounded_shortest_paths.g
    if only_walkable_road
        dist_mat = bounded_shortest_paths.dist_mat_walkable
        parents = bounded_shortest_paths.parents_walkable[from]
    else 
        dist_mat = bounded_shortest_paths.dist_mat
        parents = bounded_shortest_paths.parents[from]
    end
    if !haskey(parents, to)
        return a_star(g, from, to, dist_mat)
    end
    edges = Vector{Edge}()
    current = to 
    while current != from 
        pushfirst!(edges, Edge(parents[current], current))
        current = parents[current]
    end
    return edges
end

"""
    get_shortest_path(city_map::Map, sp_from_id, sp_to_id; only_walkable_road=false)

Return the shortest path from `from` to `to`. Has the same output as `shortest_path` from the LightOSM package 
but uses the `BoundedAllShortestPaths` cache when it exists.
"""
function get_shortest_path(city_map::Map, sp_from_id, sp_to_id; only_walkable_road=false)
    if sp_from_id == sp_to_id
        return [sp_from_id, sp_to_id]
    end
    osmgraph = city_map.graph
    graph_from_id = osmgraph.node_to_index[sp_from_id]
    graph_to_id = osmgraph.node_to_index[sp_to_id]
    
    edges = get_shortest_path(city_map.bounded_shortest_paths, graph_from_id, graph_to_id; only_walkable_road)
    # if there is no path at all
    isempty(edges) && return nothing
    osm_vertices = Vector{Int}()
    for e in edges
        push!(osm_vertices, osmgraph.index_to_node[src(e)])
    end
    push!(osm_vertices, osmgraph.index_to_node[dst(edges[end])])
    return osm_vertices
end
