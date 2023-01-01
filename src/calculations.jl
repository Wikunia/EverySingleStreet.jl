

function getxy_from_lat_lon(lat, lon, trans)
    x,y,z = trans(LLA(lat, lon))
    return x,y
end



function point_linesegment_distance(p::Point, a::Point, b::Point)
	t_hat = dot(p-a, b-a)/norm(b-a)^2
	t_star = clamp(t_hat, 0, 1)
	l(t) = a+t*(b-a)
	return norm(l(t_star)-p)
end

function get_dist_to_way(p, way::Way, trans)
    lp = Point(getxy_from_lat_lon(p.lat, p.lon, trans)...)
    min_dist = Inf
    for i in 1:length(way.nodes)-1
        w1 = way.nodes[i]
        w2 = way.nodes[i+1]
        w1p = Point(getxy_from_lat_lon(w1.lat, w1.lon, trans))
        w2p = Point(getxy_from_lat_lon(w2.lat, w2.lon, trans))
        dist = point_linesegment_distance(lp, w1p, w2p)
        min_dist = min(min_dist, dist)
    end
    return min_dist
end

function get_matching_way(map, way_ids, p)
    origin_lla = get_centroid(map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    min_dist = Inf
    min_way_id = 0 

    for way_id in way_ids
        way = map.ways[way_id]
        dist_to_way = get_dist_to_way(p, way, trans)
        if dist_to_way < min_dist
            min_way_id = way_id
            min_dist = dist_to_way
        end
    end
    if min_way_id == 0
        return nothing
    end
    return map.ways[min_way_id]
end

function get_way_kdtree(map)
    origin_lla = get_centroid(map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    n = 0
    for  way in map.ways
        for node in way.nodes
            n += 1
        end
    end
    id_to_way_id = zeros(Int, n)
    positions =  zeros(Float64, (2, n))
    i = 0
    for (way_id, way) in enumerate(map.ways)
        for node in way.nodes
            i += 1
            id_to_way_id[i] = way_id
            p = getxy_from_lat_lon(node.lat, node.lon, trans)
            positions[:, i] .= p
        end
    end
    kdtree = KDTree(positions)
    radius = longest_segment(map)*0.51
    return id_to_way_id, kdtree, radius
end

function mapping(map, paths)
    origin_lla = get_centroid(map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    id_to_way_id, way_tree, radius = get_way_kdtree(map)

    
    kps = zeros(2, )
    np = 0
    for path in paths
        np += length(path)
    end
    kps = zeros(Float64, (2, np))
    i = 0
    for path in paths
        for p in path
            i += 1
            kps[:,i] .= getxy_from_lat_lon(p.lat, p.lon, trans)
        end
    end
    @time idxs = inrange(way_tree, kps, radius)
    
    way_ids = Set()
    i = 0
    for path in paths
        for p in path
            i += 1
            search_way_ids = unique(id_to_way_id[idxs[i]])
            way = get_matching_way(map, search_way_ids, p)
            if !isnothing(way)
                push!(way_ids, way.id)
            end
        end
    end
    return way_ids
end

function longest_segment(map)
    origin_lla = get_centroid(map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    max_dist = 0.0
    for way in map.ways
        for i in 1:length(way.nodes)-1
            w1 = way.nodes[i]
            w2 = way.nodes[i+1]
            w1p = Point(getxy_from_lat_lon(w1.lat, w1.lon, trans))
            w2p = Point(getxy_from_lat_lon(w2.lat, w2.lon, trans))
            max_dist = max(max_dist, norm(w1p-w2p))
        end
    end
    return max_dist
end

function segments(way::Way)
    return [(LLA(node1.lat, node1.lon), LLA(node2.lat, node2.lon)) for (node1,node2) in zip(way.nodes[1:end-1], way.nodes[2:end])]
end

function total_length(map::Map)
    dist = 0.0
    for way in map.ways
        for segment in segments(way)
            dist += euclidean_distance(segment...)
        end
    end
    return dist
end

function total_length(paths::Vector{Vector{LLA{T}}}) where T
    dist = 0.0
    for path in paths
        for (ps, pe) in zip(path[1:end-1], path[2:end])
            dist += euclidean_distance(ps, pe)
        end
    end
    return dist
end