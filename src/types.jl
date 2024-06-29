
"""
    Node

Represents a single node in the street network.

# Fields
 - `id::Int`: Unique identifier for the node.
 - `lat::Float64`: Latitude coordinate of the node.
 - `lon::Float64``: Longitude coordinate of the node.
"""
struct Node 
    id::Int
    lat::Float64
    lon::Float64
end

Node(id, lla::LLA) = Node(id, lla.lat, lla.lon)

"""
    Way

Represents a single way (street) in the street network.
Most data is coming from openstreetmap.org so for more information about the meaning of `highway` as an example
their documentation should be checked.

# Fields
- `id::Int`: Unique identifier for the way.
- `nodes::Vector{Node}`: A vector of [`Node`](@ref) objects defining the way.
- `name::String`: Name of the way (e.g., street name).
- `highway` (String): Highway type of the way (e.g., "motorway", "residential").
- `foot` (String): Access for pedestrians (e.g., "yes", "no").
- `access` (String): General access information for the way.
- `meters::Float64`: Total length of the way in meters.
"""
struct Way
    id::Int
    nodes::Vector{Node}
    name::String
    highway::String
    foot::String
    access::String
    meters::Float64
end

"""
    HolePolygon

Represents a polygon with holes used for defining districts.

# Fields
- `outer::Vector{Point2{Float32}}`: A vector of `Point2{Float32}` objects defining the outer boundary of the polygon.
- `holes::Vector{Vector{Point2{Float32}}}`: A vector of vectors of `Point2{Float32}` objects defining any holes within the polygon.
"""
struct HolePolygon 
    outer::Vector{Point2{Float32}}
    holes::Vector{Vector{Point2{Float32}}}
end

"""
    District

Represents a district within the city.

# Fields
- `name::Symbol`: Symbolic name of the district.
- `polygons::Vector{HolePolygon}`: A vector of [`HolePolygon`](@ref) objects defining the geographical area of the district.
"""
struct District
    name::Symbol
    polygons::Vector{HolePolygon}
end

"""
    BoundedAllShortestPaths

Holds pre-computed shortest paths between nodes in a city.

# Fields
 - `g::StaticGraphs.StaticDiGraph{Int32, Int32}`: The underlying StaticDiGraph{Int32, Int32} representing the street network.
 - `dist_mat::SparseArrays.SparseMatrixCSC{Float64, Int32}`: A SparseMatrixCSC{Float64, Int32} containing the distance matrix for the nodes in a sparse matrix.
 - `distances::Vector{Dict{Int32, Float64}}`: A vector of dictionaries, where each dictionary maps a destination node ID to its shortest distance from the source node.
 - `parents::Vector{Dict{Int32, Int32}}`: A vector of dictionaries, where each dictionary maps a destination node ID to its parent node on the shortest path from the source node.
 - `distance::Float64`: The distance for which all the shortest paths are computed
"""
struct BoundedAllShortestPaths
    g::StaticGraphs.StaticDiGraph{Int32, Int32}
    dist_mat::SparseArrays.SparseMatrixCSC{Float64, Int32}
    distances::Vector{Dict{Int32, Float64}}
    parents::Vector{Dict{Int32, Int32}}
    distance::Float64
end 

"""
    AbstractSimpleMap

Abstract supertype for different map representations.
"""
abstract type AbstractSimpleMap end

"""
    NoGraphMap <: AbstractSimpleMap

Represents a map where the network graph is not explicitly stored.

# Fields
 - `osm_id_to_node_id::Dict{Int, Int}`: A dictionary mapping OpenStreetMap node IDs to internal node IDs.
 - `osm_id_to_edge_id::Dict{Int, Int}`: A dictionary mapping OpenStreetMap way (street) IDs to internal edge IDs.
 - `nodes_to_district_name::Vector{Symbol}`: A vector associating each node with the district it belongs to (represented by a symbol).
 - `nodes::Vector{Node}`: A vector of [`Node`](@ref) objects representing all nodes in the network.
 - `ways::Vector{Way}`: A vector of [`Way`](@ref) objects representing all ways (streets) in the network.
 - `walkable_road_nodes::Vector{Bool}`: A vector indicating whether each node is on a walkable road.
 - `osm_node_id_to_edge_ids::Dict{Int, Vector{Int}}`: A dictionary mapping OpenStreetMap node IDs to a vector of edge IDs it connects to.
 - `districts::Vector{District}`: A vector of [`District`](@ref) objects defining districts within the city.
"""
struct NoGraphMap <: AbstractSimpleMap
    osm_id_to_node_id::Dict{Int, Int}
    osm_id_to_edge_id::Dict{Int, Int}
    nodes_to_district_name::Vector{Symbol}
    nodes::Vector{Node}
    ways::Vector{Way}
    walkable_road_nodes::Vector{Bool}
    osm_node_id_to_edge_ids::Dict{Int, Vector{Int}}
    districts::Vector{District}
end

"""
    Map <: AbstractSimpleMap

Represents a map containing detailed street network information and potentially pre-computed shortest paths.

# Fields
 - `no_graph_map::NoGraphMap`: A [`NoGraphMap`](@ref) instance containing core street network data.
 - `graph::Union{Missing, OSMGraph})`: An optional `OSMGraph` instance representing the underlying street network graph (may be missing if shortest paths are not pre-computed).
 - `bounded_shortest_paths::Union{Missing, BoundedAllShortestPaths}`: An optional [`BoundedAllShortestPaths`](@ref) instance containing pre-computed shortest paths within a specific area (may be missing if not pre-computed).
"""
struct Map <: AbstractSimpleMap
    no_graph_map::NoGraphMap
    graph::OSMGraph
    bounded_shortest_paths::BoundedAllShortestPaths
end

"""
    GPSPoint

Represents a single point with GPS coordinates and associated time.

# Fields
 - `pos::LLA`: An LLA object containing latitude and longitude coordinates as well as altitude (not used though)
 - `time::ZonedDateTime`: The time the GPS point was recorded.
"""
struct GPSPoint 
    pos::LLA
    time::ZonedDateTime
end

"""
    Candidate

Represents a candidate point potentially lying on a walkable street segment.

# Fields
 - `measured_point::GPSPoint`: The GPS point corresponding to the candidate.
 - `lla::LLA`: The LLA object containing the candidate's latitude and longitude.
 - `way::Way`: The [`Way`](@ref) object representing the street segment the candidate is on
 - `way_is_reverse::Bool`: Indicates if the candidate is positioned along the way in the reverse direction.
 - `dist::Float64`: The distance between the measured GPS point and the candidate point on the way.
 - `λ::Float64`: The distance along the way on which the `lla` sits basically the interpolation value for the way.
"""
struct Candidate
    measured_point::GPSPoint
    lla::LLA
    way::Way
    way_is_reverse::Bool
    dist::Float64
    λ::Float64
end

"""
    GPXFile

Represents a GPX file containing a collection of GPS points.

Fields:
  - `name::String`: The name of the GPX file.
  - `gps_points::Vector{GPSPoint}`: A vector of [`GPSPoint`](@ref) objects representing the GPS points stored in the file.
"""
struct GPXFile
    name::String 
    gps_points::Vector{GPSPoint}
end

"""
    struct StreetSegment

A street segment contains of two candidates which have the following property:
They are both on the same way and have the same direction.

# Fields
- `from::Candidate` the start [`Candidate`](@ref) of the segment.
- `to::Candidate` the start [`Candidate`](@ref) of the segment.
"""
struct StreetSegment
    from::Candidate
    to::Candidate
    function StreetSegment(from, to)
        @assert from.way.id == to.way.id
        @assert from.way_is_reverse == to.way_is_reverse
        new(from, to)
    end
end

"""
    StreetPath

Represents a path through the street network defined by a sequence of connected street segments.

# Fields
  - `name::String`: A name for the path (e.g., user-defined name).
  - `subpath_id::Int`: Unique identifier for the path within a larger context (e.g., route).
  - `segments::Vector{StreetSegment}`: A vector of [`StreetSegment`](@ref) objects representing the connected street segments forming the path.
"""
struct StreetPath
    name::String
    subpath_id::Int
    segments::Vector{StreetSegment}
end

"""
    WalkedWay

Represents a way (street) segment that has been walked along, potentially with additional information about the walked parts.

Fields
  - `way::Way`: The [`Way`](@ref) object representing the underlying street segment.
  - `parts::Vector{Tuple{Float64, Float64}}`: A vector of tuples where each tuple represents a portion of the way that was walked (start and end distance along the way in meters).
"""
mutable struct WalkedWay
    way::Way
    parts::Vector{Tuple{Float64, Float64}}
end

"""
    WalkedParts

Represents a collection of walked way segments for a specific area, potentially with additional information for each way.

Fields:
  - `names::Dict{String, Vector{Int}}`: A dictionary where keys are way names (from `Way.name`) and values are vectors of integers referencing corresponding [`WalkedWay`](@ref) objects in the `ways` field.
  - `ways::Dict{Int, WalkedWay}`: A dictionary where keys are unique identifiers and values are [`WalkedWay`](@ref) objects representing the walked way segments.
"""
struct WalkedParts
    names::Dict{String, Vector{Int}} 
    ways::Dict{Int, WalkedWay}
end
