

function getxy_from_lat_lon(lat, lon, trans)
    x,y,z = trans(LLA(lat, lon))
    return x,y
end

function get_lla(p, trans)
    return trans(ENU(p...))
end


function get_interpolation_val(p::Point, a::Point, b::Point)
    t_hat = dot(p-a, b-a)/norm(b-a)^2
	t_star = clamp(t_hat, 0, 1)
    return t_star
end

get_interpolation_point(a, b, t::Float64) = a+t*(b-a)

function point_linesegment_distance(p::Point, a::Point, b::Point)
	t = get_interpolation_val(p, a, b)
    p_on_ab = get_interpolation_point(a,b, t)
	return norm(p_on_ab-p)
end

function get_candidate_on_way(p, way::Way, trans, rev_trans)
    lp = Point(getxy_from_lat_lon(p.lat, p.lon, trans)...)
    min_dist = Inf
    best_candidate = nothing
    λ = 0.0
    for i in 1:length(way.nodes)-1
        w1 = way.nodes[i]
        w2 = way.nodes[i+1]
        w1p = Point(getxy_from_lat_lon(w1.lat, w1.lon, trans))
        w2p = Point(getxy_from_lat_lon(w2.lat, w2.lon, trans))
        dist = point_linesegment_distance(lp, w1p, w2p)
        if dist < min_dist
            min_dist = dist 
            t = get_interpolation_val(lp, w1p, w2p)
            p_on_ab = get_lla(get_interpolation_point(w1p, w2p, t), rev_trans)
            best_candidate = Candidate(p_on_ab, way, dist, λ+euclidean_distance(LLA(w1.lat, w1.lon), p_on_ab))
        end
        λ += euclidean_distance(LLA(w2.lat, w2.lon), LLA(w1.lat, w1.lon))
    end
    return best_candidate
end

function get_matching_candidates(city_map, way_ids, p, origin_lla; maximum_dist=100)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)

    candidates = Vector{Candidate}()
    for way_id in way_ids
        way = city_map.ways[way_id]
        best_candidate = get_candidate_on_way(p, way, trans, rev_trans)
        if best_candidate.dist < maximum_dist
            push!(candidates, best_candidate)
        end
    end
    return candidates
end

"""
    get_way_kdtree(city_map)

Create a KDTree for all ways in the given city map.
Return 
1. a mapping from id to the way id
2. the kd tree
3. The radius for a reasonable in range search (51% of the longest line between two nodes)
"""
function get_way_kdtree(city_map)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    n = 0
    for  way in city_map.ways
        n += length(way.nodes)
    end
    id_to_way_id = zeros(Int, n)
    positions =  zeros(Float64, (2, n))
    i = 0
    for (way_id, way) in enumerate(city_map.ways)
        for node in way.nodes
            i += 1
            id_to_way_id[i] = way_id
            p = getxy_from_lat_lon(node.lat, node.lon, trans)
            positions[:, i] .= p
        end
    end
    kdtree = KDTree(positions)
    radius = longest_segment(city_map)*0.51
    return id_to_way_id, kdtree, radius
end

"""
    get_candidates(map, path)

Get a list of [`Candidate`](@ref) for each gps point in path given the underlying city_map
"""
function get_candidates(city_map, path)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    id_to_way_id, way_tree, radius = get_way_kdtree(city_map)

    kps = zeros(2, )
    np = length(path)
    kps = zeros(Float64, (2, np))
    i = 0
    for p in path
        i += 1
        kps[:,i] .= getxy_from_lat_lon(p.lat, p.lon, trans)
    end
    idxs = inrange(way_tree, kps, radius)
    
    candidates = Vector{Vector{Candidate}}()
    i = 0
    for p in path
        i += 1
        search_way_ids = unique(id_to_way_id[idxs[i]])
        p_candidates = get_matching_candidates(city_map, search_way_ids, p, origin_lla)
        push!(candidates, p_candidates)
    end
    return candidates
end

function longest_segment(city_map)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    max_dist = 0.0
    for way in city_map.ways
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

function total_length(city_map::Map)
    dist = 0.0
    for way in city_map.ways
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