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
    bounded_all_shortest_paths(osm_graph, distance, nodeid_to_local, is_walkable_road)

Call [`bounded_dijkstra`](@ref) for all nodes and return a `BoundedAllShortestPaths` object.
"""
function bounded_all_shortest_paths(osm_graph, distance, nodeid_to_local, is_walkable_road)
    g = osm_graph.graph
    dist_mat = osm_graph.weights
    all_distances = Vector{Dict{Int32, Float64}}(undef, nv(g))
    all_parents = Vector{Dict{Int32, Int32}}(undef, nv(g))
    @showprogress for i in 1:nv(g)
        if !is_walkable_road[nodeid_to_local[osm_graph.index_to_node[i]]]
            # all_distances[i] = Dict{Int32, Float64}()
            # all_parents[i] = Dict{Int32, Int32}()
            continue
        end
        distances, parents = bounded_dijkstra(g, dist_mat, i, distance)
        all_distances[i] = distances
        all_parents[i] = parents
    end
    return BoundedAllShortestPaths(g, dist_mat, all_distances, all_parents, distance)
end

"""
    get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to)

Return the same output as `a_star(g, from, to, dist_mat)` but use the cache from `BoundedAllShortestPaths`
"""
function get_shortest_path(bounded_shortest_paths::BoundedAllShortestPaths, from, to)
    g = bounded_shortest_paths.g
    dist_mat = bounded_shortest_paths.dist_mat
    parents = bounded_shortest_paths.parents[from]
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
    get_shortest_path(city_map::Map, sp_from_id, sp_to_id)

Return the shortest path from `from` to `to`. Has the same output as `shortest_path` from the LightOSM package 
but uses the `BoundedAllShortestPaths` cache when it exists.
"""
function get_shortest_path(city_map::Map, sp_from_id, sp_to_id)
    if sp_from_id == sp_to_id
        return [sp_from_id, sp_to_id]
    end
    osmgraph = city_map.graph
    graph_from_id = osmgraph.node_to_index[sp_from_id]
    graph_to_id = osmgraph.node_to_index[sp_to_id]
    
    edges = get_shortest_path(city_map.bounded_shortest_paths, graph_from_id, graph_to_id)
    # if there is no path at all
    isempty(edges) && return nothing
    osm_vertices = Vector{Int}()
    for e in edges
        push!(osm_vertices, osmgraph.index_to_node[src(e)])
    end
    push!(osm_vertices, osmgraph.index_to_node[dst(edges[end])])
    return osm_vertices
end
