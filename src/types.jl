
struct Node 
    id::Int
    lat::Float64
    lon::Float64
end

Node(id, lla::LLA) = Node(id, lla.lat, lla.lon)

struct Way
    id::Int
    nodes::Vector{Node}
    name::String
    highway::String
    foot::String
    access::String
end

struct BoundedAllShortestPaths
    g::StaticGraphs.StaticDiGraph{Int32, Int32}
    dist_mat::SparseArrays.SparseMatrixCSC{Float64, Int32}
    distances::Vector{Dict{Int32, Float64}}
    parents::Vector{Dict{Int32, Int32}}
    distance::Float64
end 

struct Map
    graph::OSMGraph
    osm_id_to_node_id::Dict{Int, Int}
    osm_id_to_edge_id::Dict{Int, Int}
    nodes::Vector{Node}
    ways::Vector{Way}
    bounded_shortest_paths::BoundedAllShortestPaths
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
