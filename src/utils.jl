
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
        JSON3.pretty(f, JSON3.write(json))
        println(f)
    end
end

function string_from_lla(lla::LLA)
    return "$(lla.lat) $(lla.lon)"
end

"""
    parse_map(fpath)

Return a [`Map`](@ref) object from the given json file path that was created using the [`download`](@ref) function.
"""
function parse_map(fpath, geojson_path=nothing)
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
    json_string = convert_keys_recursive(json)
    graph = graph_from_object(json_string; weight_type=:distance, network_type=:all)
    nodes_to_district_name = map_nodes_to_district(nodes, geojson_path)
    bounded_shortets_paths = bounded_all_shortest_paths(graph, 0.25, nodeid_to_local, walkable_road_nodes)
    return Map(graph, nodeid_to_local, wayid_to_local, nodes_to_district_name, nodes, ways, bounded_shortets_paths, walkable_road_nodes, osm_node_id_to_edge_ids)
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
    
        points = filter_path(gpxFile.tracks[1].segments[1].points, 25)
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