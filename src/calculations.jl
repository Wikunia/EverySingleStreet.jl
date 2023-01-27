

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

function get_candidate_on_way(city_map, p, way::Way, trans, rev_trans)
    lp = Point(getxy_from_lat_lon(p.lat, p.lon, trans)...)
    min_dist = Inf
    best_candidate = nothing
    cλ = 0.0
    max_λ = total_length(way)
    for i in 1:length(way.nodes)-1
        w1 = way.nodes[i]
        w2 = way.nodes[i+1]
        if !haskey(city_map.graph.nodes, w1.id) ||  !haskey(city_map.graph.nodes, w2.id)
            continue
        end
        w1p = Point(getxy_from_lat_lon(w1.lat, w1.lon, trans))
        w2p = Point(getxy_from_lat_lon(w2.lat, w2.lon, trans))
        dist = point_linesegment_distance(lp, w1p, w2p)
        if dist < min_dist
            min_dist = dist 
            t = get_interpolation_val(lp, w1p, w2p)
            p_on_ab = get_lla(get_interpolation_point(w1p, w2p, t), rev_trans)
            λ = cλ+euclidean_distance(LLA(w1.lat, w1.lon), p_on_ab)
            λ = clamp(λ, 0, max_λ)
            best_candidate = Candidate(p, p_on_ab, way, false, dist, λ)
        end
        cλ += euclidean_distance(LLA(w2.lat, w2.lon), LLA(w1.lat, w1.lon))
    end
    return best_candidate
end

function get_matching_candidates(city_map, way_ids, p, origin_lla; maximum_dist=100)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)

    candidates = Vector{Candidate}()
    for way_id in way_ids
        way = city_map.ways[way_id]
        best_candidate = get_candidate_on_way(city_map, p, way, trans, rev_trans)
        if best_candidate !== nothing && best_candidate.dist < maximum_dist
            push!(candidates, best_candidate)
            push!(candidates, get_reverse_candidate(best_candidate))
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

function filter_candidates!(candidates)
    isempty(candidates) && return candidates
    closer_dist = 25
    # if there are candidates which are closer than closer_dist0m drop all further away
    min_dist = minimum(c.dist for c in candidates)
    if min_dist > closer_dist
        return candidates
    end
    filter!(c->c.dist <= closer_dist, candidates)
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
        filter_candidates!(p_candidates)
        push!(candidates, p_candidates)
    end
    return candidates
end

"""
    get_candidate_probability(candidate::Candidate; sigma=10)

Return the emission probability of the given candidate and the standard deviation of the gps error.
"""
function get_candidate_probability(candidate::Candidate; sigma=10)
    return 1/(sqrt(2π)*sigma)*exp((-candidate.dist^2)/(2*sigma^2))
end

function get_next_node_id(candidate::Candidate)
    idx = prev_idx(candidate)+1
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return nodes[idx].id
end

function get_prev_node_id(candidate::Candidate)
    idx = prev_idx(candidate)
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return nodes[idx].id
end

function shortest_candidate_path(from::Candidate, to::Candidate, city_map)
    sp_from_id = get_next_node_id(from)
    sp_to_id = get_prev_node_id(to)
    if sp_from_id == sp_to_id
        return [sp_from_id, sp_to_id]
    end
    sp = shortest_path(city_map.graph, sp_from_id, sp_to_id)
    return sp
end

"""
    shortest_path_distance(from::Candidate, to::Candidate, city_map)

Return the shortest path distance based on equation (9) in the [fmm paper](https://www.tandfonline.com/doi/pdf/10.1080/13658816.2017.1400548?casa_token=RrxNeRXfRfkAAAAA:2IA6z4Pu-tKEcSaK44AUgQxDc-XUCBs8CSZI1qGNbUj6CpUMyA8suDUpnZ1WO3lHEUFuk1lk3s4wJtM)
"""
function shortest_path_distance(from::Candidate, to::Candidate, city_map)
    # on the same way
    if from.way.id == to.way.id
        if from.way_is_reverse != to.way_is_reverse
            return Inf
        end
        # we assume we can travel in both directions 
        d = to.λ-from.λ
        if d >= 0 
            return d
        end
    end
    # distance from end of from to begin of to
    len_way_from = total_length(from.way)
    sp = shortest_candidate_path(from, to, city_map)
    sp_from_e_to_b = total_length(city_map, sp)
    
    return  len_way_from - from.λ + sp_from_e_to_b + to.λ
end

"""
    transition_probability(from::Candidate, to::Candidate, city_map)

Return the transition_probability based on equation 10 in the  [fmm paper](https://www.tandfonline.com/doi/pdf/10.1080/13658816.2017.1400548?casa_token=RrxNeRXfRfkAAAAA:2IA6z4Pu-tKEcSaK44AUgQxDc-XUCBs8CSZI1qGNbUj6CpUMyA8suDUpnZ1WO3lHEUFuk1lk3s4wJtM)
"""
function transition_probability(from::Candidate, to::Candidate, city_map)
    d = euclidean_distance(from.measured_point, to.measured_point)
    sp_dist = shortest_path_distance(from, to, city_map)
    p = min(d, sp_dist)/max(d, sp_dist)
    if isnan(p)
        return 1.0
    end
    return p
end

function filter_path(path, dist)
    new_path = [path[1]]
    for i in 2:length(path)-1
        if euclidean_distance(new_path[end], path[i]) > dist
            push!(new_path, path[i])
        end
    end
    push!(new_path, path[end])
    saved_perc = 100*(1-(length(new_path)/length(path)))
    println("Saved $(saved_perc)%")
    return new_path
end

function get_candidates_from_idx(vec_candidates, candidate_idxs)
    candidates = Vector{Candidate}()
    ncandidates = sum(length(cands) for cands in vec_candidates)
    map_to_tpl_idx = Vector{Tuple{Int, Int}}(undef, ncandidates)
    k = 0
    for i in 1:length(vec_candidates)
        for j in 1:length(vec_candidates[i])
            k += 1
            map_to_tpl_idx[k] = (i,j)
        end
    end
    for idx in candidate_idxs
        tpl = map_to_tpl_idx[idx]
        push!(candidates, vec_candidates[tpl[1]][tpl[2]])
    end
    return candidates
end

function map_path(city_map, path)
    path = filter_path(path, 25)
    vec_candidates = EverySingleStreet.get_candidates(city_map, path)
    val, idx = findmin(length.(vec_candidates))
    if val == 0
        # Todo also do the rest later
        vec_candidates = vec_candidates[1:idx-1]
    end
    ncandidates = sum(length(cands) for cands in vec_candidates)
    @show length(vec_candidates)
    println("Average number of candidates: $(ncandidates / length(vec_candidates))")
    emission_probabilies = zeros((ncandidates, length(vec_candidates)+1))
    sidx = 1
    for (ci,candidates) in enumerate(vec_candidates)
        emission_probabilies[sidx:sidx+length(candidates)-1, ci] = EverySingleStreet.get_candidate_probability.(candidates) 
        sidx += length(candidates)
    end
    for i in 1:ncandidates
        emission_probabilies[i, end] = 1-sum(emission_probabilies[i,:])
    end
    emission_distributions = [Categorical(emission_probabilies[i,:]) for i in 1:ncandidates]
    transition_probabilities = zeros((ncandidates, ncandidates))
    sidx1 = 1
    for (ci1,candidates1) in enumerate(vec_candidates[1:end-1])
        sidx2_orig = sidx1+length(candidates1)
        candidates2 = vec_candidates[ci1+1]
        for c1 in candidates1
            sum_probs = 0.0
            sidx2 = sidx2_orig
            for c2 in candidates2
                prob = transition_probability(c1, c2, city_map)
                transition_probabilities[sidx1, sidx2] = prob
                sum_probs += prob
                sidx2 += 1
            end
            if sum_probs == 0
                transition_probabilities[sidx1,sidx2_orig:sidx2_orig+length(candidates2)-1] .= 1/length(candidates2)
            else 
                transition_probabilities[sidx1,:] ./= sum_probs
            end
            sidx1 +=1
        end
       
    end
    len_last_candidates = length(vec_candidates[end])
    for i in ncandidates-len_last_candidates+1:ncandidates
        transition_probabilities[i, i] = 1.0
    end
    println("Computed all probabilties")
    hmm = HMM(transition_probabilities, emission_distributions)
    println("Created hmm")
    best_candidates_idx = viterbi(hmm, 1:length(vec_candidates))
    best_candidates = get_candidates_from_idx(vec_candidates, best_candidates_idx)
    return best_candidates
end

function calculate_streetpath(candidates, city_map)
    segments = Vector{StreetSegment}()
    for ci in 2:length(candidates)
        current_candidate = candidates[ci-1]
        next_candidate = candidates[ci]
        if current_candidate.way.id == next_candidate.way.id
            @assert current_candidate.way_is_reverse == next_candidate.way_is_reverse
            push!(segments, StreetSegment(current_candidate, next_candidate))
        else
            sp = shortest_candidate_path(current_candidate, next_candidate, city_map)
            partial_segments = get_segments(city_map, current_candidate, next_candidate, sp)
            append!(segments, partial_segments)
        end
    end
    return StreetPath(segments)
end

function get_segments(city_map, current_candidate, next_candidate, sp)
    segments = Vector{StreetSegment}()
    if length(sp) == 2
        nodeid = get_next_node_id(current_candidate)
        node = get_node(city_map, nodeid)
        lla = LLA(node.lat, node.lon)
        c1 = Candidate(lla, lla, current_candidate.way, current_candidate.way_is_reverse, 0, total_length(current_candidate.way))
        push!(segments, StreetSegment(current_candidate, c1))
        nodeid = get_prev_node_id(next_candidate)
        node = get_node(city_map, nodeid)
        lla = LLA(node.lat, node.lon)
        c2 = Candidate(lla, lla, next_candidate.way, next_candidate.way_is_reverse, 0, total_length(next_candidate.way))
        push!(segments, StreetSegment(c2, next_candidate))
        return segments
    end

    way_segments = get_way_segments(sp, city_map)
    nodeid = get_next_node_id(current_candidate)
    node = get_node(city_map, nodeid)
    lla = LLA(node.lat, node.lon)
    c1 = Candidate(lla, lla, current_candidate.way, current_candidate.way_is_reverse, 0, total_length(current_candidate.way))
    push!(segments, StreetSegment(current_candidate, c1))

    for way_segment in way_segments
        nodes = way_segment.way.nodes
        if way_segment.rev
            nodes = reverse(nodes)
        end
        node = nodes[way_segment.from]
        lla = LLA(node.lat, node.lon)
        λ = 0.0
        if way_segment.from != 1
            λ = sum(euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1, n2) in zip(nodes[1:way_segment.from-1], nodes[2:way_segment.from]))
        end
        c1 = Candidate(lla, lla, way_segment.way, way_segment.rev, 0, λ)

        node = nodes[way_segment.to]
        lla = LLA(node.lat, node.lon)
        λ = 0.0
        if way_segment.to != 1
            λ = sum(euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1, n2) in zip(nodes[1:way_segment.to-1], nodes[2:way_segment.to]))
        end
        c2 = Candidate(lla, lla, way_segment.way, way_segment.rev, 0, λ)
        push!(segments, StreetSegment(c1, c2))
    end
    nodeid = get_prev_node_id(next_candidate)
    node = get_node(city_map, nodeid)
    lla = LLA(node.lat, node.lon)
    c2 = Candidate(lla, lla, next_candidate.way, next_candidate.way_is_reverse, 0, total_length(next_candidate.way))
    push!(segments, StreetSegment(c2, next_candidate))
    return segments
end

function get_way_segments(original_sp, city_map)
    sp = copy(original_sp)
    way_segments = Vector{NamedTuple{(:way, :rev, :from, :to), Tuple{Way, Bool, Int, Int}}}()
    while !isempty(sp)
        new_way_segment, sp = get_first_way_segment(sp, city_map)
        push!(way_segments, new_way_segment)
    end
    return way_segments
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

function total_length(city_map::Map, node_ids::Vector{Int})
    nodes = [city_map.nodes[city_map.osm_id_to_node_id[osm_id]] for osm_id in node_ids]
    dist = 0.0
    for (node1, node2) in zip(nodes[1:end-1], nodes[2:end])
        dist += euclidean_distance(LLA(node1.lat, node1.lon), LLA(node2.lat, node2.lon))
    end
    return dist
end

function total_length(way::Way)
    dist = 0.0
    for segment in segments(way)
        dist += euclidean_distance(segment...)
    end
    return dist
end

function total_length(city_map::Map)
    dist = 0.0
    for way in city_map.ways
        dist += total_length(way)
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