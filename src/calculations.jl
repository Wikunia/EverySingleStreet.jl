"""
    getxy_from_lat_lon(lat, lon, trans)

Return x, y position coordinates from lat and lon given a Geodesy trensformation.
"""
function getxy_from_lat_lon(lat, lon, trans)
    x,y,z = trans(LLA(lat, lon))
    return x,y
end

function getxy(p::GPSPoint, trans)
    x,y,z = trans(LLA(p.pos.lat, p.pos.lon))
    return x,y
end

"""
    get_lla(p, trans)

Get a LLA formatted position from a point and a Geodesy transformation.
"""
function get_lla(p, trans)
    return trans(ENU(p...))
end

"""
    get_interpolation_val(p::Point, a::Point, b::Point)

Get an interpolation value of point p between a line segment from a to b where a+t*(b-a) describes the point closes to p.
Return p which is between 0 and 1.
"""
function get_interpolation_val(p::Point, a::Point, b::Point)
    t_hat = dot(p-a, b-a)/norm(b-a)^2
	t_star = clamp(t_hat, 0, 1)
    return t_star
end

"""
    get_interpolation_point(a, b, t::Float64)

Get the value of `a+t*(b-a)`
"""
get_interpolation_point(a, b, t::Float64) = a+t*(b-a)

"""
    point_linesegment_distance(p::Point, a::Point, b::Point)

Get the smallest distance between point `p` and the line segment from `a` to `b`.
"""
function point_linesegment_distance(p::Point, a::Point, b::Point)
	t = get_interpolation_val(p, a, b)
    p_on_ab = get_interpolation_point(a,b, t)
	return norm(p_on_ab-p)
end

function get_candidate_on_way(node::Node, way::Way; rev=false)
    nodes = rev ? reverse(way.nodes) : way.nodes
    node_ids = [n.id for n in nodes]
    @assert node.id in node_ids
    lla = LLA(node.lat, node.lon)
    gpspoint = GPSPoint(lla, ZonedDateTime(now(), TimeZone("UTC")))
    λ = 0.0
    last_node = nodes[1]
    for n in nodes
        λ += euclidean_distance(LLA(last_node.lat, last_node.lon), LLA(n.lat, n.lon))
        if n.id == node.id
            break
        end
        last_node = n
    end
    return Candidate(gpspoint, lla, way, rev, 0, λ)
end

"""
    get_candidate_on_way(city_map, p, way::Way, trans, rev_trans; rev=false)

Get the best candidate of point `p` on the given `way`.
`trans` and `rev_trans` are transformations mapping from LLA to x,y and back.
`rev` can be set to true to reverse the direction of `way`.
"""
function get_candidate_on_way(city_map, p, way::Way, trans, rev_trans; rev=false)
    lp = Point(getxy(p, trans)...)
    min_dist = Inf
    best_candidate = nothing
    cλ = 0.0
    max_λ = total_length(way)
    nodes = rev ? reverse(way.nodes) : way.nodes
    for i in 1:length(nodes)-1
        w1 = nodes[i]
        w2 = nodes[i+1]
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
            best_candidate = Candidate(p, p_on_ab, way, rev, dist, λ)
        end
        cλ += euclidean_distance(LLA(w2.lat, w2.lon), LLA(w1.lat, w1.lon))
    end
    return best_candidate
end

"""
    get_matching_candidates(city_map, way_ids, p, origin_lla; maximum_dist=30)

Get matching candidates for a given point p and a list of possible way ids.
Only return candidates which have a maximum_dist (in m) to the way.
"""
function get_matching_candidates(city_map, way_ids, p, origin_lla; maximum_dist=30)
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

"""
    filter_candidates!(candidates; closer_dist=25)

Filter out candidates which are further away than `closer_dist` if there is at least candidate which is closer than `closer_dist`.
"""
function filter_candidates!(candidates; closer_dist=25)
    isempty(candidates) && return candidates
    # if there are candidates which are closer than closer_dist (in m) drop all further away
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
        kps[:,i] .= getxy(p, trans)
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
    dist = candidate.dist
    return 1/(sqrt(2π)*sigma)*exp((-dist^2)/(2*sigma^2))
end

"""
    get_next_node_id(candidate::Candidate)

Get the next node id given a candidate. This depends on `way_is_reverse` of the candidate.
Can be used to create a `StreetSegment`.
"""
function get_next_node_id(candidate::Candidate)
    idx = prev_idx(candidate)+1
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return nodes[idx].id
end

"""
    get_prev_node_id(candidate::Candidate)

Get the previous node id given a candidate. This depends on `way_is_reverse` of the candidate.
Can be used to create a `StreetSegment`.
"""
function get_prev_node_id(candidate::Candidate)
    idx = prev_idx(candidate)
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return nodes[idx].id
end

"""
    shortest_candidate_path(from::Candidate, to::Candidate, city_map)

Return the shortest path as node ids from one candidate to another given a city map.
"""
function shortest_candidate_path(from::Candidate, to::Candidate, city_map)
    sp_from_id = get_next_node_id(from)
    sp_to_id = get_prev_node_id(to)
    if sp_from_id == sp_to_id
        return [sp_from_id, sp_to_id]
    end
    sp = get_shortest_path(city_map, sp_from_id, sp_to_id)
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
    d = euclidean_distance(from.measured_point.pos, to.measured_point.pos)
    sp_dist = shortest_path_distance(from, to, city_map)
    p = min(d, sp_dist)/max(d, sp_dist)
    if isnan(p)
        return 1.0
    end
    return p
end

"""
    filter_path(gps_points, dist)

Return a new path for which the minimum distance between two consecutive points is at least `dist`.
"""
function filter_path(gps_points, dist)
    new_gps_points = [gps_points[1]]
    for i in 2:length(gps_points)-1
        if euclidean_distance(new_gps_points[end].pos, gps_points[i].pos) > dist
            push!(new_gps_points, gps_points[i])
        end
    end
    push!(new_gps_points, gps_points[end])
    saved_perc = 100*(1-(length(new_gps_points)/length(gps_points)))
    return new_gps_points
end

function filter_path(gps_points::Vector{GPX.GPXPoint}, dist)
    new_gps_points = [gps_points[1]]
    for i in 2:length(gps_points)-1
        if euclidean_distance(LLA(new_gps_points[end].lat, new_gps_points[end].lon), LLA(gps_points[i].lat, gps_points[i].lon)) > dist
            push!(new_gps_points, gps_points[i])
        end
    end
    push!(new_gps_points, gps_points[end])
    saved_perc = 100*(1-(length(new_gps_points)/length(gps_points)))
    return new_gps_points
end

"""
    get_candidates_from_idx(vec_candidates, candidate_idxs)

Return candidates from an nested vector of candiates like `[[c1,c2],[c3]]` and candiate_idxs which are 1d like `[1,3]` would return
`[c1, c3]`.
"""
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

function map_matching(city_map, gpxfile::GPXFile)
    streetpaths = Vector{StreetPath}()
    i = 1
    for (idx, candidates) in enumerate(map_path(city_map, gpxfile))
        new_streetpaths = calculate_streetpath(gpxfile.name, i, candidates, city_map)
        append!(streetpaths, new_streetpaths)
        i += length(new_streetpaths)
    end
    return streetpaths
end

"""
    map_path(city_map, gpxfile)

Map a path of gps points to the best matching candidates for the path.
1. Filter out points on the path which are closer together than 25m
2. Compute candidates for each point in the remaining path
3. Compute emission probabilties
4. Compute transition_probabilities
5. Use the Viterbi algorithm to compute the most likely path of the candidates
Return a list of lists of candidates as a path might not continously have candidates. Then it's split up into several parts.
"""
function map_path(city_map, gpxfile::GPXFile)
    gps_points = filter_path(gpxfile.gps_points, 25)
    vec_candidates = EverySingleStreet.get_candidates(city_map, gps_points)
    idxs = Vector{Tuple{Int, Int}}()
    state = :new
    current_start_idx = 0
    for (idx,candidates) in enumerate(vec_candidates)
        if state == :new
            if length(candidates) > 0
                current_start_idx = idx
                state = :ongoing
            end
        end
        if state == :ongoing
            if length(candidates) == 0
                push!(idxs, (current_start_idx, idx-1))
                current_start_idx = 0
                state = :new
            end
        end
    end
    if state == :ongoing
        push!(idxs, (current_start_idx, length(vec_candidates)))
    end
    candidates = Vector{Vector{Candidate}}()
    for tpl in idxs
        push!(candidates, map_path(city_map, vec_candidates[tpl[1]:tpl[2]]))
    end
    return candidates
end


function map_path(city_map, vec_candidates::Vector{Vector{Candidate}})
    ncandidates = sum(length(cands) for cands in vec_candidates)
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
    hmm = HMM(transition_probabilities, emission_distributions)
    best_candidates_idx = viterbi(hmm, 1:length(vec_candidates))
    best_candidates = get_candidates_from_idx(vec_candidates, best_candidates_idx)
    return best_candidates
end

"""
    calculate_streetpath(ame, subpath_id, candidates, city_map)

Generate a `StreetPath` from a list of candidates obtained by [`map_path`](@ref).
"""
function calculate_streetpath(name, subpath_id, candidates, city_map)
    streetpaths = Vector{StreetPath}()
    segments = Vector{StreetSegment}()
    for ci in 2:length(candidates)
        current_candidate = candidates[ci-1]
        next_candidate = candidates[ci]
        if current_candidate.way.id == next_candidate.way.id
            if current_candidate.way_is_reverse == next_candidate.way_is_reverse
                push!(segments, StreetSegment(current_candidate, next_candidate))
            else
                next_candidate = get_reverse_candidate(next_candidate)
                push!(segments, StreetSegment(current_candidate, next_candidate))
            end
        else
            sp = shortest_candidate_path(current_candidate, next_candidate, city_map)
            partial_segments = get_segments(city_map, current_candidate, next_candidate, sp)
            len_shortest_path =  total_length(city_map, sp)u"m"
            start_time = current_candidate.measured_point.time
            finish_time = next_candidate.measured_point.time
            duration = Quantity(finish_time - start_time)
            speed = uconvert(u"km/hr", len_shortest_path/duration)
            if speed > 20u"km/hr"
                push!(streetpaths, StreetPath(name, subpath_id, segments))
                subpath_id += 1
                segments = Vector{StreetSegment}()
            else
                append!(segments, partial_segments)
            end
        end
    end
    if !isempty(segments)
        push!(streetpaths, StreetPath(name, subpath_id, segments))
    end
    return streetpaths
end


"""
    get_segments(city_map, current_candidate, next_candidate, sp)

Given two candidates which can't be directly connected and a list of ids which form the shortest path between those candidates
this function returns a list of `StreetSegment` which connect the two candidates.
"""
function get_segments(city_map, current_candidate, next_candidate, sp)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)

    segments = Vector{StreetSegment}()
    if length(sp) == 2 && sp[1] == sp[2]
        nodeid = get_next_node_id(current_candidate)
        node = get_node(city_map, nodeid)
        gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
        c1 = get_candidate_on_way(city_map, gps_point, current_candidate.way, trans, rev_trans; rev=current_candidate.way_is_reverse)
        push!(segments, StreetSegment(current_candidate, c1))
        nodeid = get_prev_node_id(next_candidate)
        node = get_node(city_map, nodeid)
        gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
        c2 = get_candidate_on_way(city_map, gps_point, next_candidate.way, trans, rev_trans; rev=current_candidate.way_is_reverse)
        push!(segments, StreetSegment(c2, next_candidate))
        return segments
    end


    way_segments = get_way_segments(sp, city_map)
    nodeid = get_next_node_id(current_candidate)
    node = get_node(city_map, nodeid)
    gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
    c1 = get_candidate_on_way(city_map, gps_point, current_candidate.way, trans, rev_trans; rev=current_candidate.way_is_reverse)
    push!(segments, StreetSegment(current_candidate, c1))

    for way_segment in way_segments
        nodes = way_segment.way.nodes
        if way_segment.rev
            nodes = reverse(nodes)
        end
        node = nodes[way_segment.from]
        c1 = get_candidate_on_way(node, way_segment.way; rev=way_segment.rev)

        node = nodes[way_segment.to]
        c2 =  get_candidate_on_way(node, way_segment.way; rev=way_segment.rev)
        push!(segments, StreetSegment(c1, c2))
    end

    nodeid = get_prev_node_id(next_candidate)
    node = get_node(city_map, nodeid)
    gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
    c2 = get_candidate_on_way(city_map, gps_point, next_candidate.way, trans, rev_trans; rev=next_candidate.way_is_reverse)
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

function total_length(streetpath::StreetPath)
    dist = 0.0
    for segment in streetpath.segments
        dist += euclidean_distance(segment.from.measured_point.pos, segment.to.measured_point.pos)
    end
    return dist
end

function total_length(streetpaths::Vector{StreetPath})
    return sum(total_length(p) for p in streetpaths)
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

function total_length(city_map::Map; filter_fct=(way)->true)
    dist = 0.0
    for way in city_map.ways
        if filter_fct(way)
            dist += total_length(way)
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

function total_length(parts::WalkedParts; filter_fct=(way)->true)
    dist = 0.0
    for (k,way) in parts.ways
        if filter_fct(way.way)
            dist += sum(p[2]-p[1] for p in way.parts)
        end
    end
    return dist
end



function calculate_walked_parts(streetpaths::Vector{StreetPath}, city_map::Map)
    names = Dict{String, Vector{Int}}()
    ways = Dict{Int, WalkedWay}()
    for way in city_map.ways
        if haskey(names, way.name)
            push!(names[way.name], way.id)
        else
            names[way.name] = [way.id]
        end
    end

    for streetpath in streetpaths
        for segment in streetpath.segments
            len_way = total_length(segment.from.way)
            start_λ = 0.0
            finish_λ = 0.0
            if segment.from.way_is_reverse
                start_λ = len_way-segment.from.λ
                finish_λ = len_way-segment.to.λ
            else
                start_λ = segment.from.λ
                finish_λ = segment.to.λ
            end
            start_λ, finish_λ = extrema((start_λ, finish_λ))
            if start_λ == finish_λ
                continue
            end
            start_λ = clamp(start_λ-5, 0, len_way)
            finish_λ = clamp(finish_λ+5, 0, len_way)
            if !haskey(ways, segment.from.way.id)
                ways[segment.from.way.id] = WalkedWay(segment.from.way, [(start_λ, finish_λ)])
            else
                ranges = ways[segment.from.way.id].parts
                push!(ranges, (start_λ, finish_λ))
                ways = @set ways[segment.from.way.id].parts = merge_ranges(ranges)
            end
        end
    end

    return WalkedParts(names, ways)
end

function merge_ranges(ranges::Vector{Tuple{T, T}}) where T
    # sort the ranges by their start value
    sorted_ranges = sort(ranges, by = x -> x[1])

    # initialize the result with the first range
    result = [sorted_ranges[1]]

    # iterate over the rest of the ranges
    for i in 2:length(sorted_ranges)
        # if the current range overlaps with the previous range, merge them
        if sorted_ranges[i][1] <= result[end][2]
            result = @set result[end][2] = max(result[end][2], sorted_ranges[i][2])
        else
            # if there is no overlap, add the current range to the result
            push!(result, sorted_ranges[i])
        end
    end

    return result
end

function get_walked_street(walked_parts::WalkedParts, city_map::Map, name)
    way_ids = walked_parts.names[name]
    complete_dist = 0.0
    walked_dist = 0.0
    for way_id in way_ids
        local_way_id = city_map.osm_id_to_edge_id[way_id]
        way = city_map.ways[local_way_id]
        l_walked_dist, l_complete_dist = get_walked_way(walked_parts, way)
        walked_dist += l_walked_dist
        complete_dist += l_complete_dist
    end
    return walked_dist, complete_dist
end

function get_walked_way(walked_parts::WalkedParts, way)
    complete_dist = 0.0
    walked_dist = 0.0
    complete_dist += total_length(way)
    if haskey(walked_parts.ways, way.id)
        walked_way = walked_parts.ways[way.id]
        walked_dist += sum(p[2]-p[1] for p in walked_way.parts)
    end
    return walked_dist, complete_dist
end

function get_candidate_on_way(way::Way, dist)
    nodes = way.nodes
    dists = [euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1,n2) in zip(nodes[1:end-1], nodes[2:end])]
    cum_dists = cumsum(dists)
    prepend!(cum_dists, 0)
    from_idx = findlast(<=(dist), cum_dists)
    if dist >= cum_dists[end]
        lla = LLA(nodes[end].lat, nodes[end].lon)
        gpspoint = GPSPoint(lla, ZonedDateTime(now(), TimeZone("UTC")))
        candidate = Candidate(gpspoint, lla, way, false, dist, cum_dists[end])
        return candidate
    end
    to_idx = from_idx+1
    from_node = nodes[from_idx]
    to_node = nodes[to_idx]

    t = 1-(cum_dists[to_idx]-dist)/ (cum_dists[to_idx]-cum_dists[from_idx])
    lla = LLA(get_interpolation_point(Point(from_node.lat, from_node.lon), Point(to_node.lat, to_node.lon), t)...)
    gpspoint = GPSPoint(lla, ZonedDateTime(now(), TimeZone("UTC")))
    candidate = Candidate(gpspoint, lla, way, false, 0.0, dist)
    return candidate
end

function get_segments(walked_way::WalkedWay)
    segments = Vector{StreetSegment}()
    for part in walked_way.parts

        from_candidate = get_candidate_on_way(walked_way.way, part[1])
        to_candidate = get_candidate_on_way(walked_way.way, part[2])
        push!(segments, StreetSegment(from_candidate, to_candidate))
    end
    return segments
end

function get_gps_points(segment)
    prev_from_idx = prev_idx(segment.from)+1
    prev_to_idx = prev_idx(segment.to)
    points = Vector{GPSPoint}()
    if segment.from == segment.to
        return points
    end
    if prev_from_idx > prev_to_idx
        push!(points, GPSPoint(segment.from.lla , ZonedDateTime(now(), TimeZone("UTC"))))
        push!(points, GPSPoint(segment.to.lla , ZonedDateTime(now(), TimeZone("UTC"))))
        return points
    end
    points = [GPSPoint(segment.from.lla , ZonedDateTime(now(), TimeZone("UTC")))]
    nodes = segment.from.way.nodes
    if segment.from.way_is_reverse
        nodes = reverse(nodes)
    end

    for i in prev_from_idx:prev_to_idx
        node = nodes[i]
        push!(points, GPSPoint(LLA(node.lat, node.lon),  ZonedDateTime(now(), TimeZone("UTC"))))
    end
    push!(points, GPSPoint(segment.to.lla , ZonedDateTime(now(), TimeZone("UTC"))))
    return points
end

function streetpaths_to_gpx(streetpaths::Vector{StreetPath}, gps_filename)
    author = GPX.GPXAuthor("EverySingleStreet.jl")

    metadata = GPX.GPXMetadata(
        name="EverySingleStreet",
        author=author,
        time=now(localzone())
    )

    gpx = GPX.GPXDocument(metadata)

    track = GPX.new_track(gpx)

    for streetpath in streetpaths
        segments = streetpath.segments
        for segment in segments
            track_segment = GPX.new_track_segment(track)

            gps_points = get_gps_points(segment)
            for p in gps_points
                point = GPX.GPXPoint(p.pos.lat, p.pos.lon, 0, p.time, "")
                push!(track_segment, point)
            end
        end
    end

    xdoc = XMLDocument(gpx)


    save_file(xdoc, gps_filename)
    println("GPX file saved to \"$gps_filename\"")
end

function walked_parts_to_gpx(walked_parts::WalkedParts, gps_filename)
    author = GPX.GPXAuthor("EverySingleStreet.jl")

    metadata = GPX.GPXMetadata(
        name="EverySingleStreet",
        author=author,
        time=now(localzone())
    )

    gpx = GPX.GPXDocument(metadata)

    track = GPX.new_track(gpx)

    for (way_id, way) in walked_parts.ways
        segments = get_segments(way)
        for segment in segments
            track_segment = GPX.new_track_segment(track)

            gps_points = get_gps_points(segment)
            for p in gps_points
                point = GPX.GPXPoint(p.pos.lat, p.pos.lon, 0, p.time, "")
                push!(track_segment, point)
            end
        end
    end

    xdoc = XMLDocument(gpx)


    save_file(xdoc, gps_filename)
    println("GPX file saved to \"$gps_filename\"")
end

function map_matching(city_map, folder, outpath)
    streetpaths = get_gpx_paths(folder);

    save_streetpaths(outpath, Vector{StreetPath}())

    for path in streetpaths
        @show path.name
        new_streetpaths =  EverySingleStreet.map_matching(city_map, path);
        EverySingleStreet.update_streetpaths!(outpath, new_streetpaths)
    end
end

function prev_idx(nodes, λ)
    dists = [euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1,n2) in zip(nodes[1:end-1], nodes[2:end])]
    cum_dists = cumsum(dists)
    pos = findfirst(>=(λ), cum_dists)
    if pos === nothing 
        return length(cum_dists)
    end
    return pos
end

function prev_idx(candidate)
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return prev_idx(nodes, candidate.λ)
end

function get_gps_point(nodes, λ, trans, rev_trans)
    dists = [euclidean_distance(LLA(n1.lat, n1.lon), LLA(n2.lat, n2.lon)) for (n1,n2) in zip(nodes[1:end-1], nodes[2:end])]
    cum_dists = cumsum(dists)
    pos = findfirst(>=(λ), cum_dists)
    if pos === nothing 
        return LLA(nodes[end].lat, nodes[end].lon)
    end
    next_pos = pos+1
    seg_total_dist = dists[pos]
    prev_dist = pos-1 > 0 ? cum_dists[pos-1] : 0.0
    seg_dist = λ - prev_dist
    t =  seg_dist/seg_total_dist
    w1p = Point(getxy_from_lat_lon(nodes[pos].lat, nodes[pos].lon, trans))
    w2p = Point(getxy_from_lat_lon(nodes[next_pos].lat, nodes[next_pos].lon, trans))
    p_on_ab = get_lla(get_interpolation_point(w1p, w2p, t), rev_trans)

    return p_on_ab
end

function create_xml(city_map::Map, walked_parts::WalkedParts, fname)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)
    gid = 0
    ways = Vector{Way}()
    for walked_way in values(walked_parts.ways)
        nodes = walked_way.way.nodes
        for part in walked_way.parts 
            start_idx = prev_idx(nodes, part[1])
            stop_idx = prev_idx(nodes, part[2])
            part_nodes = Vector{Node}()
            gid += 1
            push!(part_nodes, Node(gid, get_gps_point(nodes, part[1], trans, rev_trans)))
            for idx in start_idx:stop_idx
                gid += 1
                push!(part_nodes, Node(gid, nodes[idx].lat, nodes[idx].lon))
            end
            gid += 1
            push!(part_nodes, Node(gid, get_gps_point(nodes, part[2], trans, rev_trans)))
            gid += 1
            push!(ways, Way(gid, part_nodes, walked_way.way.name, "walked", "yes", "yes"))
        end
    end
    
    zoned_now = ZonedDateTime(now(), TimeZone("UTC"))
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "osm")
    set_attributes(xroot, Dict("version"=>"0.6"))
    child = new_child(xroot, "bounds")
    set_attributes(child, Dict("minlon" => "7.96783", "minlat" => "53.39183", "maxlon" => "10.33196", "maxlat" => "54.06261"))

    for way in ways
        for node in way.nodes
            child = new_child(xroot, "node")
            set_attributes(child, Dict("id" => node.id, "lat" => node.lat, "lon" => node.lon, "version"=> "5", "timestamp" => zoned_now))
        end
    end
    for way in ways
        child = new_child(xroot, "way")
        set_attributes(child, Dict("id" => way.id, "version" => "1", "timestamp" => zoned_now))
        for node in way.nodes
            nd = new_child(child, "nd")
            set_attributes(nd, Dict("ref" => node.id))
        end
        tag = new_child(child, "tag")
        set_attributes(tag, Dict("k" => "name", "v" => way.name))
        tag = new_child(child, "tag")
        set_attributes(tag, Dict("k" => "highway", "v" => "primary"))
    end
    save_file(xdoc, fname)
end