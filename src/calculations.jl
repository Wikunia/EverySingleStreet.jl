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

function getxy(n::Node, trans)
    x,y,z = trans(n.lla)
    return x,y
end

function getxy(lla::LLA, trans)
    x,y,z = trans(lla)
    return x,y
end

"""
    get_lla(p, trans)

Get a LLA formatted position from a point and a Geodesy transformation.
"""
function get_lla(p, trans)
    lla = trans(ENU(p...))
    # avoid weird altitude value
    return LLA(lla.lat, lla.lon)
end

"""
    get_interpolation_val(p::Point2, a::Point2, b::Point2)

Get an interpolation value of point p between a line segment from a to b where a+t*(b-a) describes the point closes to p.
Return p which is between 0 and 1.
"""
function get_interpolation_val(p::Point2, a::Point2, b::Point2)
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
    point_linesegment_distance(p::Point2, a::Point2, b::Point2)

Get the smallest distance between point `p` and the line segment from `a` to `b`.
"""
function point_linesegment_distance(p::Point2, a::Point2, b::Point2)
	t = get_interpolation_val(p, a, b)
    p_on_ab = get_interpolation_point(a,b, t)
	return norm(p_on_ab-p)
end

"""
    node_on_way_to_candidate(idx::Int, way::Way; rev=false)
    node_on_way_to_candidate(node::Node, way::Way; rev=false)

Given either the idx of a node on a way like `way.nodes[idx]` if not reversed or the node itself
return the matching candidate by computing `λ` the time of the candidate is set to now in UTC.
"""
function node_on_way_to_candidate(idx::Int, way::Way; rev=false)
    nodes = rev ? reverse(way.nodes) : way.nodes
    node = nodes[idx]
    lla = LLA(node.lat, node.lon)
    gpspoint = GPSPoint(lla, ZonedDateTime(now(UTC), TimeZone("UTC")))
    λ = 0.0
    last_node = nodes[1]
    @views for n in nodes[1:idx]
        λ += euclidean_distance(LLA(last_node.lat, last_node.lon), LLA(n.lat, n.lon))
        last_node = n
    end
    return Candidate(gpspoint, lla, way, rev, 0, λ)
end


function node_on_way_to_candidate(node::Node, way::Way; rev=false)
    nodes = rev ? reverse(way.nodes) : way.nodes
    lla = LLA(node.lat, node.lon)
    gpspoint = GPSPoint(lla, ZonedDateTime(now(UTC), TimeZone("UTC")))
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
    lp = Point2(getxy(p, trans)...)
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
        w1p = Point2(getxy_from_lat_lon(w1.lat, w1.lon, trans))
        w2p = Point2(getxy_from_lat_lon(w2.lat, w2.lon, trans))
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
    get_matching_candidates(city_map, way_ids, p, origin_lla)

Get matching candidates for a given point p and a list of possible way ids.
Only return candidates which have a maximum_dist (in m) to the way.
"""
function get_matching_candidates(city_map, way_ids, p, origin_lla)
    maximum_dist = get_preference("CANDIDATES_MAXIMUM_DISTANCE")
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
    filter_candidates!(candidates)

Filter out candidates which are further away than `closer_dist` if there is at least candidate which is closer than `closer_dist`.
"""
function filter_candidates!(candidates)
    isempty(candidates) && return candidates
    closer_dist = get_preference("CANDIDATES_FILTER_DISTANCE")
    # if there are candidates which are closer than closer_dist (in m) drop all further away
    min_dist = minimum(c.dist for c in candidates)
    if min_dist > closer_dist
        return candidates
    end
    filter!(c->c.dist <= closer_dist, candidates)
end


"""
    find_in_between_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p1::GPSPoint, p2::GPSPoint, origin_lla::LLA, radius)

Return candidates that are as close as possible to p2 on the line from p1 to p2. 
This should be called when there are no candidates for p2 but there are for p1. 
In that case we want to find candidates for the last point that we still have candidates for.
"""
function find_in_between_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p1::GPSPoint, p2::GPSPoint, origin_lla::LLA, radius)
    dist = euclidean_distance(p1.pos, p2.pos)
    if dist < 5
        return Vector{Candidate}()
    end
    candidates = Vector{Candidate}()
    last_found_candidates = Vector{Candidate}()
    while euclidean_distance(p1.pos, p2.pos) > 5
        mid = midpoint(p1, p2)
        candidates = get_candidates(city_map, way_tree, id_to_way_id, mid, origin_lla, radius)
        if isempty(candidates)
            p2 = mid
        else 
            last_found_candidates = candidates
            p1 = mid
        end
    end
    return last_found_candidates
end


"""
    get_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p::GPSPoint, origin_lla::LLA, radius)

Get canddiates for a single point given already a kdtree which can map the point to possible candidate ways and so on.
This function is used by [`find_in_between_candidates`](@ref). 
"""
function get_candidates(city_map::AbstractSimpleMap, way_tree::KDTree, id_to_way_id, p::GPSPoint, origin_lla::LLA, radius)
    trans = ENUfromLLA(origin_lla, wgs84)
    kps = zeros(Float64, (2, 1))
    kps[:,1] .= getxy(p, trans)
    idxs = inrange(way_tree, kps, radius)[1]

    search_way_ids = unique(id_to_way_id[idxs])
    filter!(wid->iswalkable_road(city_map.ways[wid]), search_way_ids)
    p_candidates = get_matching_candidates(city_map, search_way_ids, p, origin_lla)
    filter_candidates!(p_candidates)
    return p_candidates
end

"""
    get_candidates(map, path)

Get a list of [`Candidate`](@ref) for each gps point in path given the underlying city_map
"""
function get_candidates(city_map, path)
    origin_lla = get_centroid(city_map.nodes)
    trans = ENUfromLLA(origin_lla, wgs84)
    id_to_way_id, way_tree, radius = get_way_kdtree(city_map)
    radius += get_preference("CANDIDATES_MAXIMUM_DISTANCE")

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
    previously_found = false
    for p in path
        i += 1
        search_way_ids = unique(id_to_way_id[idxs[i]])
        filter!(wid->iswalkable_road(city_map.ways[wid]), search_way_ids)
        p_candidates = get_matching_candidates(city_map, search_way_ids, p, origin_lla)
        filter_candidates!(p_candidates)
        if isempty(p_candidates)
            if previously_found
                previously_found = false
                push!(candidates, find_in_between_candidates(city_map, way_tree, id_to_way_id, path[i-1], p, origin_lla, radius))
            end
        else
            if !previously_found && i > 1 # find a point closest to the previous which didn't have any canddiates
                push!(candidates, find_in_between_candidates(city_map, way_tree, id_to_way_id, p, path[i-1], origin_lla, radius))
            end
            previously_found = true
            push!(candidates, p_candidates)
        end
    end
    return candidates
end

"""
    get_candidate_probability(candidate::Candidate)

Return the emission probability of the given candidate and the standard deviation of the gps error.
"""
function get_candidate_probability(candidate::Candidate)
    dist = candidate.dist
    sigma = get_preference("GPS_STD_DEV")
    return 1/(sqrt(2π)*sigma)*exp((-dist^2)/(2*sigma^2))
end

"""
    get_next_node_id(candidate::Candidate)

Get the next node id given a candidate. This depends on `way_is_reverse` of the candidate.
Can be used to create a `StreetSegment`.
"""
function get_next_node_id(candidate::Candidate)
    @assert iswalkable_road(candidate.way)
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
    @assert iswalkable_road(candidate.way)
    idx = prev_idx(candidate)
    nodes = candidate.way.nodes
    if candidate.way_is_reverse
        nodes = reverse(nodes)
    end
    return nodes[idx].id
end

"""
    shortest_candidate_path(from::Candidate, to::Candidate, city_map; only_walkable_road=false)

Return the shortest path as node ids from one candidate to another given a city map.
"""
function shortest_candidate_path(from::Candidate, to::Candidate, city_map; only_walkable_road=false)
    from_ids = (get_prev_node_id(from), get_next_node_id(from))
    to_ids = (get_prev_node_id(to), get_next_node_id(to))
    shortest_dist = Inf
    shortest_sp = nothing
    for sp_from_id in from_ids
        for sp_to_id in to_ids
            sp = get_shortest_path(city_map, sp_from_id, sp_to_id; only_walkable_road)
            isnothing(sp) && continue
            dist = total_length(city_map, sp)
            dist += euclidean_distance(from.lla, get_lla(city_map, sp[1]))
            dist += euclidean_distance(get_lla(city_map, sp[end]), to.lla)
            if dist < shortest_dist 
                shortest_dist = dist
                shortest_sp = sp
            end
        end
    end
    return shortest_sp
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
    isnothing(sp) && return Inf
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

function create_local_json(city_map, node_ids_start, fpath)
    way_ids = Vector{Int}()
    osm_node_ids = [city_map.nodes[i].id for i in node_ids_start]
    for osm_node_id in osm_node_ids
        append!(way_ids, city_map.osm_node_id_to_edge_ids[osm_node_id])
    end
    unique!(way_ids)
    ways = [city_map.ways[way_id] for way_id in way_ids]
    osm_node_ids = Set{Int}()
    for way in ways
        union!(osm_node_ids, [n.id for n in way.nodes])
    end

    nodes = [city_map.nodes[city_map.osm_id_to_node_id[i]] for i in osm_node_ids]
    json = Dict{Symbol, Any}()
    json[:version] = 0.6
    json[:elements] = Vector{Dict{Symbol, Any}}()
    element = Dict{Symbol, Any}()
    element[:id] = 0
    element[:type] = "count"
    element[:tags] = Dict{Symbol, Any}()
    element[:tags][:relations] = "0"
    element[:tags][:nodes] = "$(length(nodes))"
    element[:tags][:ways]  = "$(length(ways))"
    element[:tags][:total] = "$(length(nodes)+length(ways))"
    push!(json[:elements], element)

    for node in nodes 
        element = Dict{Symbol, Any}()
        element[:id] = node.id
        element[:type] = "node"
        element[:lat] = node.lat
        element[:lon] = node.lon
        push!(json[:elements], element)
    end

    for way in ways 
        element = Dict{Symbol, Any}()
        element[:id] = way.id
        element[:type] = "way"
        element[:nodes] = [node.id for node in way.nodes]
        element[:tags] = Dict{Symbol, String}()
        element[:tags][:name] = way.name
        element[:tags][:highway] = way.highway
        element[:tags][:foot] = way.foot
        element[:tags][:access] = way.access

        push!(json[:elements], element)
    end

    open(fpath, "w") do f
        JSON3.write(f, json)
        println(f)
    end
end

function get_local_map(city_map, gps_points, map_local_path)
    padding=get_preference("LOCAL_MAP_PADDING")
   
    origin_lla = get_centroid(city_map.nodes)
    kd_tree = get_kd_tree_from_points(gps_points, origin_lla)

    node_transformed_points = get_transformed_points(city_map.nodes, origin_lla)
    _, dists = nn(kd_tree, node_transformed_points)
    node_ids = findall(<=(padding), dists)
    create_local_json(city_map, node_ids, map_local_path)

    return parse_map(map_local_path)
end

function map_matching(fpath, city_map::AbstractSimpleMap, walked_parts::WalkedParts, map_local_path="tmp_local_map.json")
    name = basename(fpath)
    gps_points = get_gps_points(fpath)
    return map_matching(gps_points, city_map, walked_parts, map_local_path; name)
end

function map_matching(gps_points::Vector{GPSPoint}, city_map::AbstractSimpleMap, walked_parts::WalkedParts, map_local_path="tmp_local_map.json"; name=nothing)
    if isnothing(name)
        name = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    end
    map_local = get_local_map(city_map, gps_points, map_local_path)
    if isempty(map_local.ways)
        return (this_walked_parts = WalkedParts(), walked_parts = walked_parts, added_kms = 0.0, this_walked_road_km = 0.0)     
    end
    streetpaths = map_matching(map_local, name, gps_points)
    prev_walked_road_km = total_length(walked_parts; filter_fct=(way)->EverySingleStreet.iswalkable_road(way))/1000
    this_walked_parts = calculate_walked_parts(streetpaths, city_map)
    walked_parts = calculate_walked_parts(streetpaths, city_map, walked_parts.ways)
    now_walked_road_km = total_length(walked_parts; filter_fct=(way)->EverySingleStreet.iswalkable_road(way))/1000
    this_walked_road_km = total_length(this_walked_parts; filter_fct=(way)->EverySingleStreet.iswalkable_road(way))/1000
    return (walked_parts = walked_parts, this_walked_parts = this_walked_parts, added_kms = now_walked_road_km - prev_walked_road_km, this_walked_road_km = this_walked_road_km)
end

function map_matching(city_map, gpxfile::GPXFile)
    return map_matching(city_map, gpxfile.name, gpxfile.gps_points)
end

function map_matching(city_map, gpx_name, gpxpoints::Vector{GPSPoint}; point_dist=5, allow_recursive=true)
    streetpaths = Vector{StreetPath}()
    i = 1
    candidates_vec = map_path(city_map, gpxpoints; point_dist)
    for (idx, candidates) in enumerate(candidates_vec)
        new_streetpaths = calculate_streetpath(gpx_name, i, candidates, city_map; allow_recursive)
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
    return map_path(city_map, gpxfile.gps_points)
end

function map_path(city_map, gps_points::Vector{GPSPoint}; point_dist=5)  
    gps_points = simplify(gps_points, point_dist)
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
    calculate_streetpath(name, subpath_id, candidates, city_map)

Generate a `StreetPath` from a list of candidates obtained by [`map_path`](@ref).
"""
function calculate_streetpath(name, subpath_id, candidates, city_map; allow_recursive=true)
    streetpaths = Vector{StreetPath}()
    segments = Vector{StreetSegment}()
    for ci in 2:length(candidates)
        current_candidate = candidates[ci-1]
        next_candidate = candidates[ci]
        
        sp = shortest_candidate_path(current_candidate, next_candidate, city_map)
        if isnothing(sp)
            if !isempty(segments)
                push!(streetpaths, StreetPath(name, subpath_id, segments))
            subpath_id += 1
                segments = Vector{StreetSegment}()
            end
            continue
        end
        len_shortest_path =  total_length(city_map, sp)u"m"
        any_non_walkable_road = uses_any_non_walkable_road(city_map, sp)
        start_time = current_candidate.measured_point.time
        finish_time = next_candidate.measured_point.time
        duration = Quantity(finish_time - start_time)

        speed = uconvert(u"km/hr", len_shortest_path/duration)

        # don't check speed limit if the added amount is short
        too_fast = speed > 20u"km/hr" && len_shortest_path > 50u"m"

        if !too_fast && any_non_walkable_road
            sp_only_walkable = shortest_candidate_path(current_candidate, next_candidate, city_map; only_walkable_road=true)
            if !isnothing(sp_only_walkable)
                len_sp_only_walkable =  total_length(city_map, sp_only_walkable)u"m"
                speed_only_walkable = uconvert(u"km/hr", len_sp_only_walkable/duration)
                too_fast_only_walkable = speed_only_walkable > 20u"km/hr" && len_sp_only_walkable > 50u"m"
                if !too_fast_only_walkable && len_sp_only_walkable - 2*get_preference("EXTEND_WALKED_WAY_UP_TO")u"m" <= len_shortest_path 
                    sp = sp_only_walkable
                    any_non_walkable_road = false
                end
            end
        end
        partial_segments = get_segments(city_map, current_candidate, next_candidate, sp)
        
        # if recursive is allowed and some of the shortest path are not via walkable roads 
        # try to match those again to walkable roads
        # also check if the shortest path is already too fast. If that is the case then this doesn't need to be computed
        if !too_fast && allow_recursive && any_non_walkable_road
            if !isempty(segments)
                push!(streetpaths, StreetPath(name, subpath_id, segments))
                subpath_id += 1
                segments = Vector{StreetSegment}()
            end
            nodes = [city_map.nodes[city_map.osm_id_to_node_id[nid]] for nid in sp]
            points = [LLA(n.lat, n.lon) for n in nodes]
            pushfirst!(points, current_candidate.measured_point.pos)
            push!(points, next_candidate.measured_point.pos)
            start_time = current_candidate.measured_point.time
            finish_time = next_candidate.measured_point.time
            duration = Quantity(finish_time - start_time)
            len_shortest_path =  total_length(city_map, sp)u"m"
            distances = [euclidean_distance(p1, p2)u"m" for (p1, p2) in zip(points[1:end-1], points[2:end])]
            cum_dists = cumsum(distances)
            pushfirst!(cum_dists, 0.0u"m")
            times = start_time .+ ceil.(Second, [d/len_shortest_path * duration for d in cum_dists])
            gps_points = [GPSPoint(p, t) for (p, t) in zip(points, times)]
            extra_paths = map_matching(city_map, "extra_$ci", gps_points; allow_recursive=false, point_dist=0.0)
            for extra_path in extra_paths
                push!(streetpaths, StreetPath(name, subpath_id, extra_path.segments))
                subpath_id += 1
            end
            continue
        end
        
        if too_fast && !isempty(segments)
            push!(streetpaths, StreetPath(name, subpath_id, segments))
            subpath_id += 1
            segments = Vector{StreetSegment}()
        else
            append!(segments, partial_segments)
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
        nodeid = sp[1]
        node = get_node(city_map, nodeid)
        gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
        c1 = get_candidate_on_way(city_map, gps_point, current_candidate.way, trans, rev_trans; rev=current_candidate.way_is_reverse)
        push!(segments, StreetSegment(current_candidate, c1))
        node = get_node(city_map, nodeid)
        gps_point = GPSPoint(LLA(node.lat, node.lon), current_candidate.measured_point.time)
        c2 = get_candidate_on_way(city_map, gps_point, next_candidate.way, trans, rev_trans; rev=next_candidate.way_is_reverse)
        push!(segments, StreetSegment(c2, next_candidate))
        return segments
    end


    way_segments = get_way_segments(sp, city_map)
    nodeid = sp[1]
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
        c1 = node_on_way_to_candidate(node, way_segment.way; rev=way_segment.rev)

        node = nodes[way_segment.to]
        c2 =  node_on_way_to_candidate(node, way_segment.way; rev=way_segment.rev)
        push!(segments, StreetSegment(c1, c2))
    end

    nodeid = sp[end]
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
            w1p = Point2(getxy_from_lat_lon(w1.lat, w1.lon, trans))
            w2p = Point2(getxy_from_lat_lon(w2.lat, w2.lon, trans))
            max_dist = max(max_dist, norm(w1p-w2p))
        end
    end
    return max_dist
end

function segments(way::Way)
    return segments(way.nodes)
end


function segments(nodes::Vector{Node})
    return [(LLA(node1.lat, node1.lon), LLA(node2.lat, node2.lon)) for (node1,node2) in zip(nodes[1:end-1], nodes[2:end])]
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
    return way.meters
end

function total_length(nodes)
    dist = 0.0
    for segment in segments(nodes)
        dist += euclidean_distance(segment...)
    end
    return dist
end

function total_length(city_map::AbstractSimpleMap; filter_fct=(way)->true)
    return total_length(city_map.ways; filter_fct)
end

function total_length(ways::Vector{Way}; filter_fct=(way)->true)
    dist = 0.0
    for way in ways
        if filter_fct(way)
            dist += total_length(way)
        end
    end
    return dist
end

function total_length(path::Vector{LLA{T}}) where T
    dist = 0.0
    for (ps, pe) in zip(path[1:end-1], path[2:end])
        dist += euclidean_distance(ps, pe)
    end
    return dist
end

function total_length(paths::Vector{Vector{LLA{T}}}) where T
    return sum(total_length(path for path in paths))
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

function streetpaths_to_walked_parts(streetpaths::Vector{StreetPath}, city_ways::Vector{Way}, walked_ways=Dict{Int, WalkedWay}())
    names = Dict{String, Vector{Int}}()
    for way in city_ways
        if haskey(names, way.name)
            push!(names[way.name], way.id)
        else
            names[way.name] = [way.id]
        end
    end
    possible_way_ids = Set(way.id for way in city_ways if iswalkable_road(way))

    for streetpath in streetpaths
        for segment in streetpath.segments
            segment.from.way.id in possible_way_ids || continue
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
            if !haskey(walked_ways, segment.from.way.id)
                walked_ways[segment.from.way.id] = WalkedWay(segment.from.way, [(start_λ, finish_λ)])
            else
                ranges = walked_ways[segment.from.way.id].parts
                push!(ranges, (start_λ, finish_λ))
                walked_ways[segment.from.way.id].parts = merge_ranges(ranges)
            end
        end
    end
    return WalkedParts(names, walked_ways)
end

"""
    calculate_walked_parts(streetpaths::Vector{StreetPath}, city_map::AbstractSimpleMap, walked_ways=Dict{Int, WalkedWay}())

Return `WalkedParts` given the streetpath segments and the possible ways of the city.
Can be added to already existing `walked_ways`. Filters everything which isn't a `walkable_road` inside `city_ways`.   
Calls several function to get a better approximation of what was actually walked like 
    - add extra buffer at start and end of streets to not miss the last few meters of a dead end street as an example [`extend_walked_parts!`](@ref)
    - closes circles like roundabouts, parts at ends of some dead end streets
"""
function calculate_walked_parts(streetpaths::Vector{StreetPath}, city_map::AbstractSimpleMap, walked_ways=Dict{Int, WalkedWay}())
    walked_parts = streetpaths_to_walked_parts(streetpaths, city_map.ways, walked_ways)
    extend_walked_parts!(walked_parts, city_map)
    return walked_parts
end

function extend_walked_parts!(walked_parts, city_map)
    extend_walked_parts_connectors!(walked_parts, city_map)
    extend_walked_parts_simple!(walked_parts)
    extend_walked_parts_cycle!(walked_parts)
end

"""
    extend_walked_parts_connectors!(walked_parts, city_map)

1. Find all the city ways that are shorter than 2*EXTEND_WALKED_WAY_UP_TO
2. Check if the end of those ways is part of some already walked 
"""
function extend_walked_parts_connectors!(walked_parts, city_map)
    max_connector_length = 2*get_preference("EXTEND_WALKED_WAY_UP_TO")
    city_ways = city_map.ways
    walkable_roads = filter(iswalkable_road, city_ways)
    connectors = filter(w->w.meters < max_connector_length, walkable_roads)
    filter!(w->!isfullywalked(walked_parts,w), connectors)
    node_tpls = [(w.nodes[1], w.nodes[end]) for w in connectors]
    added_m = 0.0
    for (connector, node_tpl) in zip(connectors, node_tpls)
        if check_if_visited_node(city_map, node_tpl[1], walked_parts.ways)
            if check_if_visited_node(city_map, node_tpl[2], walked_parts.ways)
                if haskey(walked_parts.ways, connector.id)
                    added_m -= sum(p->p[2]-p[1], walked_parts.ways[connector.id].parts) 
                end
                walked_parts.ways[connector.id] = WalkedWay(connector, [(0.0, connector.meters)])
                added_m += connector.meters
            end
        end
    end
end

"""
    isfullywalked(walked_parts, way::Way)

Return whether the way is already completely walked.
"""
function isfullywalked(walked_parts, way::Way)
    ww = get(walked_parts.ways, way.id, nothing)
    isnothing(ww) && return false 
    length(ww.parts) > 1 && return false 
    part = ww.parts[1]
    return part[1] == 0.0 && part[2] ≈ way.meters
end

"""
    check_if_visited_node(city_map, node::Node, walked_parts)

Return whether the node was visited already.
"""
function check_if_visited_node(city_map, node::Node, walked_parts)
    pws = get_possible_ways(city_map, node.id)
    for pw in pws 
        c = node_on_way_to_candidate(node, pw)
        ww = get(walked_parts, c.way.id, nothing)
        isnothing(ww) && continue
        any(p->p[1] <= c.λ <= p[2], ww.parts) && return true
    end
    return false
end


"""
    extend_walked_parts_simple!(walked_parts)

Extend the walked parts by adding up to `EXTEND_WALKED_WAY_UP_TO` to the start and end of a walked way.
As well as to fill gaps in between parts.
"""
function extend_walked_parts_simple!(walked_parts)
    extend_up_to = get_preference("EXTEND_WALKED_WAY_UP_TO")
    for (idx, wway) in walked_parts.ways
        len_way = total_length(wway.way)
        new_parts = Vector{Tuple{Float64, Float64}}()
        # extend gaps 
        for (part1, part2) in zip(wway.parts[1:end-1], wway.parts[2:end])
            s,e = part1[1], part1[2]
            if part1[2] + extend_up_to >= part2[1]
                e = part2[1]
            end
            push!(new_parts, (s,e))
        end
        push!(new_parts, wway.parts[end])
        extended_parts = copy(new_parts)
        new_parts = Vector{Tuple{Float64, Float64}}()
        for part in extended_parts
            s = part[1]
            e = part[2]
            if part[1] < extend_up_to
                s = 0.0
            end
            if part[2] > len_way - extend_up_to
                e = len_way
            end
            push!(new_parts, (s,e))
        end
        new_parts = merge_ranges(new_parts)
        walked_parts.ways[idx].parts = new_parts
    end
    return walked_parts
end

"""
    extend_walked_parts_cycle!(walked_parts)

Extend the walked parts by finishing a cycle if more than `MIN_FILL_CYCLE_PERC` perc of a cycle of length maximum `MAX_FILL_CYCLE_LENGTH` is walked.
"""
function extend_walked_parts_cycle!(walked_parts)
    max_cycle_len = get_preference("MAX_FILL_CYCLE_LENGTH")
    min_fill_perc = get_preference("MIN_FILL_CYCLE_PERC")
    for (idx, wway) in walked_parts.ways
        hascycle(wway.way) || continue
        # get cycles as parts of way
        way = wway.way 
        g = get_directed_graph(way)
        cycles = simplecycles(g)
        cycles_λ = Vector{Tuple{Float64, Float64}}()
        for cycle in cycles
            c1 = node_on_way_to_candidate(cycle[1], way)
            c2 = node_on_way_to_candidate(cycle[end], way)
            if c2.λ - c1.λ < max_cycle_len
                push!(cycles_λ, (c1.λ, c2.λ))
            end
        end
        new_parts = copy(wway.parts)
        for part in wway.parts 
            for cycle_λ in cycles_λ
                s = max(cycle_λ[1], part[1]) 
                e = min(cycle_λ[2], part[2]) 
                if s < e 
                    cycle_part_len = e - s 
                    cycle_len = cycle_λ[2] - cycle_λ[1]
                    if cycle_part_len/cycle_len > min_fill_perc/100
                        push!(new_parts, cycle_λ)  
                    end
                end
            end
        end
        new_parts = merge_ranges(new_parts)
        walked_parts.ways[idx].parts = new_parts
    end
    return walked_parts
end

"""
    intersection(range1::Tuple{T, T}, range2::Tuple{T, T}) where T

Return the intersection of the two ranges as a `Tuple{T, T}` or `nothing` if they don't intersect.
"""
function intersection(range1::Tuple{T, T}, range2::Tuple{T, T}) where T
    intersection_start = max(range1[1], range2[1])
    intersection_end = min(range1[2], range2[2])

    if intersection_start <= intersection_end  # Check for actual intersection
       return (intersection_start, intersection_end)
    else 
        return nothing
    end
end

"""
    merge_ranges(ranges::Vector{Tuple{T, T}}) where T

Merge the given ranges together and return a new vector of merged ranges which are also sorted.
"""
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

function get_segments(walked_way::WalkedWay)
    segments = Vector{StreetSegment}()
    for part in walked_way.parts

        from_candidate = node_on_way_to_candidate(walked_way.way, part[1])
        to_candidate = node_on_way_to_candidate(walked_way.way, part[2])
        push!(segments, StreetSegment(from_candidate, to_candidate))
    end
    return segments
end

function get_gps_points(segment::StreetSegment)
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
    LightXML.free(xdoc)
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
    LightXML.free(xdoc)
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

function update_map_matching(city_map, folder, jld2_path)
    streetpaths = get_gpx_paths(folder);
    saved_streetpaths = Vector{StreetPath}()
    if isfile(jld2_path)
        saved_streetpaths = load(jld2_path)["streetpaths"]
    end
    saved_path_names = Set(streetpath.name for streetpath in saved_streetpaths)
    # @show sort(collect(saved_path_names))
    for path in streetpaths
        path.name in saved_path_names && continue
        @show path.name
        new_streetpaths =  EverySingleStreet.map_matching(city_map, path);
        append!(saved_streetpaths, new_streetpaths)
        save_streetpaths(jld2_path, saved_streetpaths)
    end
    return saved_streetpaths
end

function full_update_routine(city_map, folder, jld2_path, output_folder=".")
    streetpaths = update_map_matching(city_map, folder, jld2_path)
    walked_parts = calculate_walked_parts(streetpaths, city_map);
    xml_path = joinpath(output_folder, "walked.xml")
    osm_path = joinpath(output_folder, "walked.osm.pbf")
    create_xml(city_map.nodes, walked_parts, xml_path)
    run(`osmosis --read-xml $xml_path --write-pbf $osm_path`)
    return walked_parts
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
    w1p = Point2(getxy_from_lat_lon(nodes[pos].lat, nodes[pos].lon, trans))
    w2p = Point2(getxy_from_lat_lon(nodes[next_pos].lat, nodes[next_pos].lon, trans))
    p_on_ab = get_lla(get_interpolation_point(w1p, w2p, t), rev_trans)

    return p_on_ab
end

function get_lla(city_map::AbstractSimpleMap, osm_node_id::Int)
    node = city_map.nodes[city_map.osm_id_to_node_id[osm_node_id]]
    return LLA(node.lat, node.lon)
end

function get_ways_from_walked_parts(walked_parts::WalkedParts)
    origin_lla = get_centroid(walked_parts)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)
    gid = 0
    ways = Vector{Way}()
    for walked_way in values(walked_parts.ways)
        nodes = walked_way.way.nodes
        for part in walked_way.parts 
            start_idx = prev_idx(nodes, part[1])+1
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
            push!(ways, Way(gid, part_nodes, walked_way.way.name, "walked", "yes", "yes", total_length(part_nodes)))
        end
    end
    return ways
end

function create_xml(all_nodes::Vector{Node}, walked_parts::WalkedParts, fname; districts=Vector{District}(), district_levels=Dict{Symbol, Int}())
    ways = get_ways_from_walked_parts(walked_parts)
    gid = ways[end].id

    minlat, maxlat = extrema(n.lat for n in all_nodes) 
    minlon, maxlon = extrema(n.lon for n in all_nodes) 
    
    zoned_now = ZonedDateTime(now(), TimeZone("UTC"))
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "osm")
    set_attributes(xroot, Dict("version"=>"0.6"))
    child = new_child(xroot, "bounds")
    set_attributes(child, Dict("minlon" => minlon, "minlat" => minlat, "maxlon" => maxlon, "maxlat" => maxlat))

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

    for district in districts
        for hpolygon in district.polygons
            start_gid = gid+1
            for pos in hpolygon.outer
                gid += 1
                child = new_child(xroot, "node")
                set_attributes(child, Dict("id" => gid, "lat" => pos[2], "lon" => pos[1], "version"=> "5", "timestamp" => zoned_now))
            end
            ids_outer = start_gid:gid
            ids_holes = Vector{UnitRange{Int}}()
            for hole in hpolygon.holes
                start_gid = gid+1
                for pos in hole
                    gid += 1
                    child = new_child(xroot, "node")
                    set_attributes(child, Dict("id" => gid, "lat" => pos[2], "lon" => pos[1], "version"=> "5", "timestamp" => zoned_now))
                end
                push!(ids_holes, start_gid:gid)
            end
            
            gid += 1
            relation = new_child(xroot, "relation")
            set_attributes(relation, Dict("id" => gid, "version"=> "5", "timestamp" => zoned_now))
            tag = new_child(relation, "tag")
            set_attributes(tag, Dict("k" => "type", "v" => "boundary"))
            
            tag = new_child(relation, "tag")
            set_attributes(tag, Dict("k" => "admin_level", "v" => "9"))

            tag = new_child(relation, "tag")
            set_attributes(tag, Dict("k" => "boundary", "v" => "administrative"))

            if haskey(district_levels, district.name)
                for district_level in district_levels[district.name] 
                    tag = new_child(relation, "tag")
                    set_attributes(tag, Dict("k" => "name", "v" => district_level))
                end
            end
            
            member = new_child(relation, "member")
            set_attributes(member, Dict("type" => "way", "ref" => gid+1, "role" => "outer"))
            for hi in 1:length(hpolygon.holes)
                member = new_child(relation, "member")
                set_attributes(member, Dict("type" => "way", "ref" => gid+1+hi, "role" => "inner"))
            end

            way = new_child(xroot, "way")
            set_attributes(way, Dict("type" => "way", "id" => gid+1, "version"=> "5", "timestamp" => zoned_now))
            for id in ids_outer
                nd = new_child(way, "nd")
                set_attributes(nd, Dict("ref" => id))
            end
            if !isempty(ids_outer)
                nd = new_child(way, "nd")
                set_attributes(nd, Dict("ref" => ids_outer[1]))
            end

            for hi in  1:length(hpolygon.holes)
                way = new_child(xroot, "way")
                set_attributes(way, Dict("type" => "way", "id" => gid+1+hi, "version"=> "5", "timestamp" => zoned_now))
                for id in ids_holes[hi]
                    nd = new_child(way, "nd")
                    set_attributes(nd, Dict("ref" => id))
                end
                if !isempty(ids_holes[hi])
                    nd = new_child(way, "nd")
                    set_attributes(nd, Dict("ref" => ids_holes[hi][1]))
                end
            end
        end
    end

    save_file(xdoc, fname)
    LightXML.free(xdoc)
end