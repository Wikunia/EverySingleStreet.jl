
"""
    download(place_name, filepath)  

Download the road network of the given place and write it as a json file to the given filepath
"""
function download(place_name, filepath)
    download_osm_network(
        :place_name;
        network_type=:drive,
        place_name,
        save_to_file_location=filepath
    )
end

function string_from_lla(lla::LLA)
    return "$(lla.lat) $(lla.lon)"
end

"""
    parse_map(fpath)

Return a [`Map`](@ref) object from the given json file path that was created using the [`download`](@ref) function.
"""
function parse_map(fpath)
    json = readjson(fpath)
    elements = json[:elements]
    counter = elements[1]
    @assert counter[:type] == "count"
    nodes = Vector{Node}(undef, parse(Int, counter[:tags][:nodes]))
    ways = Vector{Way}(undef, parse(Int, counter[:tags][:ways]))
    nodeid_to_local = Dict{Int, Int}()
    wayid_to_local = Dict{Int, Int}()
    node_counter = 1
    way_counter = 1
    for element in elements
        if element[:type] == "node"
            nodeid_to_local[element[:id]] = node_counter
            nodes[node_counter] = Node(element[:id], element[:lat], element[:lon])
            node_counter += 1
        end
        if element[:type] == "way"
            element[:tags][:oneway] = "no"
            way_nodes = [nodes[nodeid_to_local[node_id]] for node_id in element[:nodes]]
            wayid_to_local[element[:id]] = way_counter
            ways[way_counter] = Way(element[:id], way_nodes, get(element[:tags],:name, ""), get(element[:tags],:highway, ""), get(element[:tags],:foot, ""), get(element[:tags],:access, ""))
            way_counter += 1
        end
    end
    json_string = convert_keys_recursive(json)
    graph = graph_from_object(json_string; weight_type=:distance, network_type=:all)
    return Map(graph, nodeid_to_local, wayid_to_local, nodes, ways)
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

function get_first_way_segment(sp, city_map::Map)
    best_way_segment = nothing
    best_len = 0
    for (rev,func) in zip([false, true], [identity, reverse])
        for way in city_map.ways
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