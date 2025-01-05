
"""
    download(place_name, filepath)  

Download the road network of the given place and write it as a json file to the given filepath
"""
function download(place_name, filepath)
    download_osm_network(
        :place_name;
        network_type=:all,
        place_name,
        save_to_file_location=filepath
    )
end

function download(south_west::LLA, north_east::LLA, filepath)
    minlat = south_west.lat
    minlon = south_west.lon
    maxlat = north_east.lat
    maxlon = north_east.lon
    download_osm_network(
        :bbox;
        network_type=:all,
        minlat, minlon, maxlat, maxlon,
        save_to_file_location=filepath
    )
end

function filter_walkable_json!(filepath)
    json = readjson(filepath)
    new_elements = Vector{Dict}()
    removed_nways = 0
    for element in json[:elements]
        if element[:type] != "way"
            push!(new_elements, element)
            continue
        end
        new_way = Way(
            element[:id], Vector{Node}(), 
            get(element[:tags],:name, ""), 
            get(element[:tags],:highway, ""), 
            get(element[:tags],:foot, ""), 
            get(element[:tags],:access, ""),
            0.0
        )
        if iswalkable(new_way)
            push!(new_elements, element)
            continue 
        end
        removed_nways += 1
    end
    json[:elements] = new_elements

    open(filepath, "w") do f
        JSON3.write(f, json)
        println(f)
    end
end

function string_from_lla(lla::LLA)
    return "$(lla.lat) $(lla.lon)"
end


function parse_no_graph_map(fpath, geojson_path=nothing)
    json = readjson(fpath)
    elements = json[:elements]
    counter = elements[1]
    @assert counter[:type] == "count"
    nodes = Vector{Node}(undef, parse(Int, counter[:tags][:nodes]))
    ways = Vector{Way}()
    nodeid_to_local = Dict{Int, Int}()
    wayid_to_local = Dict{Int, Int}()
    osm_node_id_to_edge_ids = Dict{Int, Vector{Int}}()
    node_counter = 1
    way_counter = 1
    encountered_first_way = false
    walkable_road_nodes = [false]
    for element in elements
        if element[:type] == "node"
            @assert !encountered_first_way
            nodeid_to_local[element[:id]] = node_counter
            nodes[node_counter] = Node(element[:id], element[:lat], element[:lon])
            osm_node_id_to_edge_ids[element[:id]] = Vector{Int}()
            node_counter += 1
        end
        if element[:type] == "way"
            if !encountered_first_way
                walkable_road_nodes = zeros(Bool, node_counter-1)
                encountered_first_way = true
            end
            element[:tags][:oneway] = "no"
            way_nodes = [nodes[nodeid_to_local[node_id]] for node_id in element[:nodes]]
            local_node_ids = [nodeid_to_local[node_id] for node_id in element[:nodes]]
            new_way = Way(
                element[:id], way_nodes, 
                get(element[:tags],:name, ""), 
                get(element[:tags],:highway, ""), 
                get(element[:tags],:foot, ""), 
                get(element[:tags],:access, ""),
                total_length(way_nodes)
            )
            for node_id in element[:nodes]
                push!(osm_node_id_to_edge_ids[node_id], way_counter)
            end
            wayid_to_local[element[:id]] = way_counter
            if iswalkable_road(new_way)
                walkable_road_nodes[local_node_ids] .= true
            end
            push!(ways, new_way)
            way_counter += 1
        end
    end
    
    districts = get_districts(geojson_path) 
    nodes_to_district_name = map_nodes_to_district(nodes, geojson_path)
    return json, NoGraphMap(nodeid_to_local, wayid_to_local, nodes_to_district_name, nodes, ways, walkable_road_nodes, osm_node_id_to_edge_ids, districts)
end

"""
    parse_map(fpath)

Return a [`Map`](@ref) object from the given json file path that was created using the [`download`](@ref) function.
"""
function parse_map(fpath, geojson_path=nothing)
    json, no_graph_map = parse_no_graph_map(fpath, geojson_path)
    json_string = convert_keys_recursive(json)
    graph = graph_from_object(json_string; weight_type=:distance, network_type=:all, largest_connected_component=false)
    bounded_shortest_paths = bounded_all_shortest_paths(graph, 0.25, no_graph_map.osm_id_to_node_id, no_graph_map.walkable_road_nodes)
    return Map(no_graph_map, graph, bounded_shortest_paths)
end

function parse_gpx(fpath)
    gpxFile = GPX.read_gpx_file(fpath)
    @assert length(gpxFile.tracks) == 1
    @assert length(gpxFile.tracks[1].segments) == 1

    points = [GPSPoint(LLA(p.lat, p.lon, p.ele), p.time) for p in gpxFile.tracks[1].segments[1].points]
    filename = splitext(basename(fpath))[1]
    return GPXFile(filename, points)
end


function combine_gpx_tracks(folder)
    author = GPX.GPXAuthor("EverySingleStreet.jl")

    metadata = GPX.GPXMetadata(
        name="EverySingleStreet",
        author=author,
        time=now(localzone())
    )

    gpx = GPX.GPXDocument(metadata)

    track = GPX.new_track(gpx)
    


    for fname in readdir(folder)
        splitext(fname)[2] != ".gpx" && continue
        track_segment = GPX.new_track_segment(track)

        gpxFile = GPX.read_gpx_file(joinpath(folder,fname))
    
        @assert length(gpxFile.tracks) == 1
        @assert length(gpxFile.tracks[1].segments) == 1
    
        points = simplify(gpxFile.tracks[1].segments[1].points, 5)
        for p in points
            point = GPX.GPXPoint(p.lat, p.lon, p.ele, p.time, p.desc)
            push!(track_segment, point)
        end
    end

    xdoc = XMLDocument(gpx)


    fname = "generated.gpx"
    save_file(xdoc, fname)
    println("GPX file saved to \"$fname\"")
end

function readjson(fpath)
    json_string = read(fpath, String)
    return copy(JSON3.read(json_string))
end

function get_gpx_paths(folder)
    return [parse_gpx(joinpath(folder,fname)) for fname in readdir(folder) if splitext(fname)[2] == ".gpx"]
end

function get_centroid(walked_parts::WalkedParts)
    node_ids = Set{Int}()
    nodes = Vector{Node}()
    for (name, walked_way) in walked_parts.ways
        for node in walked_way.way.nodes 
            node.id in node_ids && continue
            push!(node_ids, node.id)
            push!(nodes, node)
        end
    end 
    return get_centroid(nodes)
end

function get_centroid(points::Vector{GPSPoint})
    lon = mean(point.pos.lon for point in points)
    lat = mean(point.pos.lat for point in points)
    return LLA(lat, lon)
end

function get_centroid(nodes::Vector{Node})
    lon = mean(node.lon for node in nodes)
    lat = mean(node.lat for node in nodes)
    return LLA(lat, lon)
end

function convert_keys_recursive(d::Dict)
    new_d = Dict{String,Any}()
    for (k,v) in d
        new_d[string(k)] = convert_keys_recursive(v)
    end
    return new_d
end

function convert_keys_recursive(v::Vector)
    return [convert_keys_recursive(d) for d in v]
end

convert_keys_recursive(x) = x

function get_reverse_candidate(candidate::Candidate)
    max_len = total_length(candidate.way)
    λ = clamp(max_len-candidate.λ, 0, max_len)
    return Candidate(
        candidate.measured_point, 
        candidate.lla,
        candidate.way,
        !candidate.way_is_reverse,
        candidate.dist,
        λ
    )
end

function get_node(city_map::Map, nodeid)
    return city_map.nodes[city_map.osm_id_to_node_id[nodeid]]
end

function get_possible_ways(city_map, node_id)
    way_ids = city_map.osm_node_id_to_edge_ids[node_id]
    return city_map.ways[way_ids]
end

"""
    get_first_way_segment(sp, city_map::Map)

Return a way that includes the longest first part of the given shortest path 
as well as the remaining shortest path which isn't part of this way.
The return format consists of:
- A named tuple describing the best way for the first segment: `(way::Way, rev::Bool, from::Int, to::Int)`
- The remaining shortest path nodes
"""
function get_first_way_segment(sp, city_map::Map)
    best_way_segment = nothing
    possible_ways = get_possible_ways(city_map, sp[1])
    best_len = 0
    for (rev,func) in zip([false, true], [identity, reverse])
        for way in possible_ways
            node_ids = func([n.id for n in way.nodes])
            last_idx = 0
            while last_idx !== nothing
                start_pos_idx = findnext(==(sp[1]), node_ids, last_idx+1)
                pos_idx = start_pos_idx
                spi = 1
                if pos_idx !== nothing && pos_idx != length(node_ids)
                    while pos_idx+1 <= length(node_ids) && spi+1 <= length(sp) && node_ids[pos_idx+1] == sp[spi+1]
                        len = spi+1
                        pos_idx += 1
                        spi += 1
                        if len > best_len
                            best_len = len
                            best_way_segment = (way=way, rev=rev, from=start_pos_idx, to=pos_idx)
                        end
                    end
                end
                last_idx = pos_idx
            end
        end
    end
    if !isnothing(best_way_segment)
        # by default always have the last already found still in it to not have a gap
        new_sp = sp[best_len:end]
        # return either length at least two or an empty one
        if length(new_sp) == 1
            return best_way_segment, Int[]
        end
        return best_way_segment, new_sp
    end
    error("Couldn't find which contains at least the first two points of $sp directly after another")
end

function save_streetpaths(filename, streetpaths::Vector{StreetPath})
    save(filename, Dict("streetpaths" => streetpaths))
end

function update_streetpaths!(filename, update_streetpaths::Vector{StreetPath}; verbose=false)
    streetpaths = load(filename, "streetpaths")
    for streetpath in update_streetpaths
        replaced = false
        for (idx,saved_streetpath) in enumerate(streetpaths)
            if streetpath.name == saved_streetpath.name && streetpath.subpath_id == saved_streetpath.subpath_id
                streetpaths[idx] = streetpath
                verbose && println("Replaced streetpaths with name $(streetpath.name)")
                replaced = true
                break
            end
        end
        replaced && continue
        verbose && println("Added streetpaths with name $(streetpath.name)")
        push!(streetpaths, streetpath)
    end
    save_streetpaths(filename, streetpaths)
end


function add_streetpaths!(filename, added_streetpaths::Vector{StreetPath})
    streetpaths = load(filename, "streetpaths")
    append!(streetpaths, added_streetpaths)
    save_streetpaths(filename, streetpaths)
end

function iswalkable_road(way::Way)
    way.highway in ["", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "living_street", "pedestrian", "secondary_link"] || return false
    way.foot in ["no", "private", "discouraged"] && return false
    way.access in ["no", "private", "customers"] && return false 
    return true 
end

function iswalkable(way::Way)
    iswalkable_road(way) && return true
    return way.highway in ["footway", "track", "steps", "path", "service"]
end

function get_gps_points(fpath)
    strava_json = readjson(fpath)
    start_time = ZonedDateTime(DateTime(strava_json[:start_time][1:end-1], "yyyy-mm-ddTHH:MM:SS"), tz"UTC")
    gpx_points = Vector{GPSPoint}()
    for (tdelta, latlon) in zip(strava_json[:times], strava_json[:latlon])
        t = start_time+Second(tdelta)
        gpx_point = GPSPoint(LLA(latlon[1], latlon[2]), t)
        push!(gpx_points, gpx_point)
    end
    return gpx_points
end

function bbox(points::Vector{GPSPoint}, padding_m=0)
    minlat, maxlat = extrema(p.pos.lat for p in points) 
    minlon, maxlon = extrema(p.pos.lon for p in points) 
    origin_lla = LLA((maxlat+minlat)/2, (minlon+maxlon)/2)
    south_west = LLA(minlat, minlon)
    north_east = LLA(maxlat, maxlon)
    trans = ENUfromLLA(origin_lla, wgs84)
    rev_trans = LLAfromENU(origin_lla, wgs84)
    south_west = get_lla(trans(south_west)[1:2] .- [padding_m, padding_m], rev_trans)
    north_east = get_lla(trans(north_east)[1:2] .+ [padding_m, padding_m], rev_trans)
    return (south_west = south_west, north_east = north_east)
end

@inline function Base.getproperty(map::Map, s::Symbol)
    if s in (
        :osm_id_to_node_id,
        :osm_id_to_edge_id,
        :nodes_to_district_name,
        :nodes,
        :ways,
        :walkable_road_nodes,
        :osm_node_id_to_edge_ids
    )
        Core.getproperty(Core.getproperty(map, :no_graph_map), s)
    else
        getfield(map, s)
    end
end

"""
Use a non-recursive Douglas-Peucker algorithm to simplify a polygon. Used by `simplify()`.

    douglas_peucker(pointlist::Array, start_index, last_index, epsilon)
"""
function douglas_peucker(pointlist::Vector, start_index, last_index, epsilon)
    temp_stack = Tuple{Int, Int}[]
    push!(temp_stack, (start_index, last_index))
    global_start_index = start_index
    keep_list = trues(length(pointlist))
    while length(temp_stack) > 0
        start_index = first(temp_stack[end])
        last_index =  last(temp_stack[end])
        pop!(temp_stack)
        dmax = 0.0
        index = start_index
        for i in index + 1:last_index - 1
            if (keep_list[i - global_start_index])
                d = pointlinedistance(pointlist[i], pointlist[start_index], pointlist[last_index])
                if d > dmax
                    index = i
                    dmax = d
                end
            end
        end
        if dmax > epsilon
            push!(temp_stack, (start_index, index))
            push!(temp_stack, (index, last_index))
        else
            for i in start_index + 2:last_index - 1 # 2 seems to keep the starting point...
                keep_list[i - global_start_index] = false
            end
        end
    end
    return pointlist[keep_list]
end

"""
Simplify a polygon:

    simplify(pointlist::Array, detail=0.1)

`detail` is the maximum approximation error of simplified polygon.
"""
function simplify(pointlist::Vector{GPSPoint}, detail=0.1)
    isempty(pointlist) && return pointlist
    origin_lla = get_centroid(pointlist)
    trans = ENUfromLLA(origin_lla, wgs84)
    transformed_points = Vector{Point2{Float64}}()
    for point in pointlist
        transformed_point = Point2(getxy_from_lat_lon(point.pos.lat, point.pos.lon, trans))
        push!(transformed_points, transformed_point)
    end

    simplified = douglas_peucker(transformed_points, 1, length(transformed_points), detail)
    gps_points = Vector{GPSPoint}()
    for p in simplified 
        idx = findfirst(==(p), transformed_points)
        @assert !isnothing(idx)
        push!(gps_points, pointlist[idx])
    end
    return gps_points
end

"""
    pointlinedistance(p::Point2, a::Point2, b::Point2)

Find the distance between a point `p` and a line between two points `a` and `b`.
"""
function pointlinedistance(p::Point2, a::Point2, b::Point2)
  dx = b[1] - a[1]
  dy = b[2] - a[2]
  return abs(p[1] * dy - p[2] * dx + b[1] * a[2] - b[2] * a[1]) / hypot(dx, dy);
end

function set_preferences!(pairs...)
    @set_preferences!(pairs...)
end

function get_preference(key)
    return @load_preference(key, DEFAULT_PREFERNECE_VALUE[key])
end
                                                                                                                                                                                                
hascycle(way::Way) = !allunique(way.nodes)

function get_directed_graph(way::Way)
    g = SimpleDiGraph(length(way.nodes))
    for i in 1:length(way.nodes)-1
        add_edge!(g, i, i+1)
    end
    for (nidx,node) in enumerate(way.nodes)
        idxs = findall(n->(n.id == node.id), way.nodes)
        filter!(i->i < nidx , idxs)
        for idx in idxs
            add_edge!(g, nidx, idx)
        end
    end
    return g
end
   
"""
    points2geojson(vec_candidates::Vector{Vector{Candidate}}, geojson_path)
    points2geojson(points::Vector{<:LLA}, geojson_path)

Debug functions to visualize candidates or simple LLA points using tools like: geojson.io/
"""
function points2geojson(vec_candidates::Vector{Vector{Candidate}}, geojson_path)
    llas = LLA[]
    for vec_candidate in vec_candidates
        append!(llas, [c.lla for c in vec_candidate])
    end
    return points2geojson(llas, geojson_path)
end

function points2geojson(points::Vector{<:LLA}, geojson_path)
    # Create GeoJSON features from points
    features = [
        Dict(
            "type" => "Feature",
            "geometry" => Dict(
                "type" => "Point",
                "coordinates" => [point.lon, point.lat]
            ),
            "properties" => Dict(
                "id" => i 
            )
        ) for (i,point) in enumerate(points)
    ]

    # Create a FeatureCollection
    geojson = Dict(
        "type" => "FeatureCollection",
        "features" => features
    )

    # Write to a GeoJSON file
    open(geojson_path, "w") do file
        JSON3.write(file, geojson)
    end
end

function ways2geojson(ways, geojson_path)
    # Create GeoJSON features from points
    features = [
        Dict(
            "type" => "Feature",
            "geometry" => Dict(
                "type" => "LineString",
                "coordinates" => [[node.lon, node.lat] for node in way.nodes]
            ),
            "properties" => Dict(
                "id" => i 
            )
        ) for (i,way) in enumerate(ways)
    ]

    # Create a FeatureCollection
    geojson = Dict(
        "type" => "FeatureCollection",
        "features" => features
    )

    # Write to a GeoJSON file
    open(geojson_path, "w") do file
        JSON3.write(file, geojson)
    end
end

"""
    midpoint(p1::GPSPoint, p2::GPSPoint)

Calculates the midpoint between two GPS points.

This function takes two `GPSPoint` objects as input, each containing a position
represented by `LLA` (Latitude, Longitude, Altitude) coordinates and a
timestamp represented by `ZonedDateTime`. It calculates the midpoint in
ECEF (Earth-Centered, Earth-Fixed) coordinates for accuracy, then converts
the result back to LLA. The time of the midpoint is also interpolated linearly.

# Arguments
- `p1::GPSPoint`: The first GPS point.
- `p2::GPSPoint`: The second GPS point.

# Returns
A `GPSPoint` object representing the midpoint between `p1` and `p2`.
"""
function midpoint(p1::GPSPoint, p2::GPSPoint)
    ecef1 = ECEF(p1.pos, wgs84)
    ecef2 = ECEF(p2.pos, wgs84)

    mid_ecef = (ecef1 + ecef2) / 2

    mid_lla = LLA(mid_ecef, wgs84)

    mid_time = p1.time + (p2.time - p1.time) ÷ 2

    return GPSPoint(mid_lla, mid_time)
end