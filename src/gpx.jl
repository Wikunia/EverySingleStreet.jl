"""
Copied from: https://github.com/scls19fr/GPX.jl/blob/master/src/GPX.jl (MIT License)
Module for working with [GPX file format](https://en.wikipedia.org/wiki/GPS_Exchange_Format).

Be aware that this module don't fully follow 1.1 version of this format.
"""
module GPX

using LightXML
using Dates
using TimeZones

import Base: push!, iterate, length, getindex

export GPXAuthor, GPXMetadata, GPXPoint, 
    GPXTrackSegment, GPXDocument, GPXRoute
export read_gpx_file, parse_gpx_string
export new_track, new_track_segment, last_segment

# not exported but very useful: GPXTracks

"""
Constant defining GPX Format version
"""
const GPX_VERSION = "1.1"

"""
Constant dictionary defining default namespaces
"""
const GPX_NS = Dict(
    "" => "http://www.topografix.com/GPX/1/1"
)
# http://www.topografix.com/GPX/1/1/gpx.xsd

"""
Constant defining default GPX file creator name
"""
const GPX_CREATOR = "GPX.jl"


"""
    GPXAuthor(name; email="", link="")

The person or organization who created the GPX file.
"""
struct GPXAuthor
    name::String
    email::String
	link::String
end
function GPXAuthor(name; email="", link="")
    GPXAuthor(name, email, link)
end

"""
    GPXMetadata(; name="", author=GPXAuthor(""), time=now(tz"UTC"))

Information about the GPX file, author, and copyright restrictions goes in the metadata section.
Providing rich, meaningful information about your GPX files allows others to search for and use your GPS data.

# Parameters
- `name`: The name of the GPX file.
- `desc`: A description of the contents of the GPX file.
- `author`: The person or organization who created the GPX file.
- `copyright`: Copyright and license information governing use of the file.
- `link`: URLs associated with the location described in the file.
- `time`: The creation date of the file.
- `keywords`: Keywords associated with the file. Search engines or databases can use this information to classify the data.
- `bounds`: Minimum and maximum coordinates which describe the extent of the coordinates in the file.
- `extensions`:
"""
struct GPXMetadata
    name::String
    desc::String
    author::GPXAuthor
    copyright::String
    link::String
    time::ZonedDateTime
    keywords::String
    bounds::String
end
GPXMetadata(; name="", desc="", author=GPXAuthor(""), copyright="", link="", time=now(tz"UTC"), keywords="", bounds="") = GPXMetadata(name, desc, author, copyright, link, time, keywords, bounds)


"""
    GPXPoint(lat::Float64, lon::Float64, ele::Float64, time::ZonedDateTime, desc::String)

A geographic point with optional elevation and time. Available for use by other schemas.

A Track Point holds the coordinates, elevation, timestamp, and metadata for a single point in a track.

# Parameters
- `lat`: The latitude of the point. This is always in decimal degrees, and always in WGS84 datum.
- `lon`: The longitude of the point. This is always in decimal degrees, and always in WGS84 datum.
- `ele`: The elevation (in meters) of the point.
- `time`: The time that the point was recorded.
- `desc`: A text description of the element.
"""
struct GPXPoint
    lat::Float64
    lon::Float64
    ele::Float64
    time::ZonedDateTime
    desc::String
end

"""
    GPXTrackSegment()

A Track Segment holds a list of Track Points which are logically connected in order. To represent a single GPS track where GPS reception was lost, or the GPS receiver was turned off, start a new Track Segment for each continuous span of track data
"""
struct GPXTrackSegment
    points::Vector{GPXPoint}
    GPXTrackSegment() = new(GPXPoint[])
end

"""
    GPXTrack()

A Track holds a list of Track Segments.
"""
struct GPXTrack
    segments::Vector{GPXTrackSegment}
    GPXTrack() = new(GPXTrackSegment[])
end

"""
    GPXTracks()

A list of tracks.
"""
struct GPXTracks
    collection::Vector{GPXTrack}
    GPXTracks() = new(Vector{GPXTrack}[])
end
push!(tracks::GPXTracks, track::GPXTrack) = push!(tracks.collection, track)
iterate(tracks::GPXTracks) = iterate(tracks.collection)
iterate(tracks::GPXTracks, i) = iterate(tracks.collection, i)
length(tracks::GPXTracks) = length(tracks.collection)
getindex(tracks::GPXTracks, i) = getindex(tracks.collection, i)

"""

ToDo
"""
struct GPXRoute
end


"""
    GPXDocument(metadata::GPXMetadata; tracks=GPXTracks(), namespaces=GPX_NS, version=GPX_VERSION, creator=GPX_CREATOR)

GPX documents contain a metadata header, followed by waypoints, routes, and tracks.
You can add your own elements to the extensions section of the GPX document.    
"""
struct GPXDocument
    namespaces::Dict{String, String}

    version::String
    creator::String
    metadata::GPXMetadata

	waypoints::Vector{GPXPoint}
	routes::Vector{GPXRoute}
    tracks::GPXTracks

    function GPXDocument(metadata::GPXMetadata; tracks=GPXTracks(), namespaces=GPX_NS, version=GPX_VERSION, creator=GPX_CREATOR)
        waypoints = GPXPoint[]
        routes = GPXRoute[]
        new(namespaces, version, creator, metadata, waypoints, routes, tracks)
    end
end

"""
    new_track(gpx::GPXDocument) -> GPXTrack

Get a new track from a `GPXDocument`.  
"""
function new_track(gpx::GPXDocument)
    return new_track(gpx.tracks)
end

"""
    new_track(tracks::GPXTracks) -> GPXTrack

Get a new track from a collection of tracks (ie `GPXTracks`).
"""
function new_track(tracks::GPXTracks)
    track = GPXTrack()
    push!(tracks, track)
    return track
end

"""
    new_track_segment(track::GPXTrack) -> GPXTrackSegment

Get a new track segment from a `GPXTrack`.
"""
function new_track_segment(track::GPXTrack)
    track_segment = GPXTrackSegment()
    push!(track.segments, track_segment)
    return track_segment
end

"""
    push!(track_segment::GPXTrackSegment, point::GPXPoint)

Insert one `GPXPoint` at the end of a `GPXTrackSegment`.
"""
function push!(track_segment::GPXTrackSegment, point::GPXPoint)
    push!(track_segment.points, point)
end

"""
    last_segment(gpx) -> GPXTrackSegment
    

"""
function last_segment(gpx::GPXDocument)
    track = isempty(gpx.tracks) ? new_track(gpx) : gpx.tracks[end]
    segment = isempty(track.segments) ? new_track_segment(track) : track.segments[end]
    return segment
end


"""
    XMLDocument(gpx::GPXDocument)

Create a `LightXML.XMLDocument` from a `GPXDocument`.

`XMLDocument` can then be saved to file or printed.
"""
function LightXML.XMLDocument(gpx::GPXDocument)
    xdoc = LightXML.XMLDocument()

    xroot = create_root(xdoc, "gpx")
    for (ns_key, ns_url) in gpx.namespaces
        key = ns_key == "" ? "xmlns" : "xmlns:$ns_key"
        set_attribute(xroot, key, ns_url)
    end
    set_attribute(xroot, "version", gpx.version)
    set_attribute(xroot, "creator", gpx.creator)
    
    x_metadata = new_child(xroot, "metadata")
    x_metadata_name = new_child(x_metadata, "name")
    add_text(x_metadata_name, gpx.metadata.name)
    x_metadata_author = new_child(x_metadata, "author")
    x_metadata_author_name = new_child(x_metadata_author, "name")
    add_text(x_metadata_author_name, gpx.metadata.author.name)
    
    for track in gpx.tracks
        x_trk = new_child(xroot, "trk")
        for segment in track.segments
            x_trkseg = new_child(x_trk, "trkseg")
            for point in segment.points
                x_trkpt = new_child(x_trkseg, "trkpt")
                lat, lon = point.lat, point.lon
                set_attribute(x_trkpt, "lat", "$lat")
                set_attribute(x_trkpt, "lon", "$lon")
                x_ele = new_child(x_trkpt, "ele")
                ele = point.ele
                add_text(x_ele, "$ele")
                x_time = new_child(x_trkpt, "time")
                dt = point.time
                s_dt_utc = "$(dt.utc_datetime)Z"
                add_text(x_time, s_dt_utc)
                desc = point.desc
                if desc != ""
                    x_desc = new_child(x_trkpt, "desc")
                    add_text(x_desc, desc)
                end
            end
        end
    end
    
    return xdoc
end

"""
    read_gpx_file(fname) -> GPXDocument

Read GPX file from filename `fname`.
"""
function read_gpx_file(fname)
    xdoc = parse_file(fname)
    return _parse_gpx(xdoc)
end

"""
    parse_gpx_string(s) -> GPXDocument

Parse GPX data from String `s`.
"""
function parse_gpx_string(s)
    xdoc = parse_string(s)
    return _parse_gpx(xdoc)
end


"""
    _parse_gpx(xdoc::XMLDocument) -> GPXDocument

Parse `XMLDocument` and return a `GPXDocument`.
"""
function _parse_gpx(xdoc::XMLDocument)
    gpxs = root(xdoc)

    metadata = GPXMetadata(
        name="title",
        author=GPXAuthor("author"),
        time=ZonedDateTime("2019-01-01T00:00:00.000+00:00"),
    )

    tracks = GPXTracks()
    
    for x_gpx in child_elements(gpxs)
        if name(x_gpx) == "trk"
            track = new_track(tracks)
            for x_track in child_elements(x_gpx)
                if name(x_track) == "trkseg"
                    track_segment = new_track_segment(track)
                    for x_track_segment in child_elements(x_track)
                        if name(x_track_segment) == "trkpt"
                            d_att = attributes_dict(x_track_segment)
                            dt = ZonedDateTime(0, tz"UTC")
                            lat = parse(Float64, d_att["lat"])
                            lon = parse(Float64, d_att["lon"])
                            ele = 0.0
                            desc = ""
                            for x_track_point in child_elements(x_track_segment)
                                if name(x_track_point) == "time"
                                    s = content(x_track_point)
                                    # Required for optional millisecond field
                                    # see https://github.com/JuliaTime/TimeZones.jl/issues/83
                                    # Check for period ('.') as indication of
                                    # datetime format
                                    dt = if occursin('.', s)
                                        parse(ZonedDateTime, s, dateformat"yyyy-mm-ddTHH:MM:SS.ssszzz")
                                    else
                                        parse(ZonedDateTime, s, dateformat"yyyy-mm-ddTHH:MM:SSzzz")
                                    end
                                elseif name(x_track_point) == "ele"
                                    s = content(x_track_point)
                                    ele = parse(Float64, s)
                                elseif name(x_track_point) == "desc"
                                    desc = content(x_track_point)
                                    # desc = ""  # for debug
                                end
                            end
                            point = GPXPoint(lat, lon, ele, dt, desc)
                            push!(track_segment, point)
                        end
                    end
                end
            end
        end
    end

    free(xdoc)

    gpx = GPXDocument(metadata, tracks=tracks)

    return gpx
end

end # module