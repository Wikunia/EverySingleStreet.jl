
struct Node 
    id::Int
    lat::Float64
    lon::Float64
end

struct Way
    id::Int
    nodes::Vector{Node}
    name::String
    highway::String
    foot::String
    access::String
end

struct HolePolygon 
    outer::Vector{Point2{Float32}}
    holes::Vector{Vector{Point2{Float32}}}
end

struct District
    name::Symbol
    polygons::Vector{HolePolygon}
end

struct Map
    graph::OSMGraph
    osm_id_to_node_id::Dict{Int, Int}
    osm_id_to_edge_id::Dict{Int, Int}
    nodes_to_district_name::Vector{Symbol}
    nodes::Vector{Node}
    ways::Vector{Way}
end

struct GPSPoint 
    pos::LLA
    time::ZonedDateTime
end
struct Candidate
    measured_point::GPSPoint
    lla::LLA
    way::Way
    way_is_reverse::Bool
    dist::Float64
    λ::Float64
end

struct GPXFile
    name::String 
    gps_points::Vector{GPSPoint}
end

"""
    struct StreetSegment

A street segment contains of two candidates which have the following property:
They are both on the same way and have the same direction.
"""
struct StreetSegment
    from::Candidate
    to::Candidate
end

struct StreetPath
    name::String
    subpath_id::Int
    segments::Vector{StreetSegment}
end

struct WalkedWay
    way::Way
    parts::Vector{Tuple{Float64, Float64}}
end

struct WalkedParts
    names::Dict{String, Vector{Int}} 
    ways::Dict{Int, WalkedWay}
end