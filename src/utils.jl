
"""
    download(place_name, filepath)  

Download the road network of the given place and write it as a json file to the given filepath
"""
function download(place_name, filepath)
    download_osm_network(
        :place_name;
        place_name,
        save_to_file_location=filepath
    )
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
            way_nodes = [nodes[nodeid_to_local[node_id]] for node_id in element[:nodes]]
            wayid_to_local[element[:id]] = way_counter
            ways[way_counter] = Way(element[:id], way_nodes, get(element[:tags],:name, ""))
            way_counter += 1
        end
    end
    return Map(nodes, ways)
end

function parse_gpx(fpath)
    gpxFile = read_gpx_file(fpath)
    @assert length(gpxFile.tracks) == 1
    @assert length(gpxFile.tracks[1].segments) == 1

    return [LLA(p.lat, p.lon, p.ele) for p in gpxFile.tracks[1].segments[1].points]
end


function combine_gpx_tracks(folder)
    author = GPXAuthor("EverySingleStreet.jl")

    metadata = GPXMetadata(
        name="07/11/2019 LFBI (09:32) LFBI (11:34)",
        author=author,
        time=now(localzone())
    )

    gpx = GPXDocument(metadata)

    track = new_track(gpx)
    


    for fname in readdir(folder)
        splitext(fname)[2] != ".gpx" && continue
        track_segment = new_track_segment(track)

        gpxFile = read_gpx_file(joinpath(folder,fname))
        @assert length(gpxFile.tracks) == 1
        @assert length(gpxFile.tracks[1].segments) == 1
    
        for p in gpxFile.tracks[1].segments[1].points
            point = GPXPoint(p.lat, p.lon, p.ele, p.time, p.desc)
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