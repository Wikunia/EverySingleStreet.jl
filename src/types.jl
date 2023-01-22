
struct Node 
    id::Int
    lat::Float64
    lon::Float64
end

struct Way
    id::Int
    nodes::Vector{Node}
    name::String
end

struct Map
    graph::OSMGraph
    osm_id_to_node_id::Dict{Int, Int}
    osm_id_to_edge_id::Dict{Int, Int}
    nodes::Vector{Node}
    ways::Vector{Way}
end

struct Candidate
    measured_point::LLA
    lla::LLA
    way::Way
    way_is_reverse::Bool
    dist::Float64
    λ::Float64
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
    segments::Vector{StreetSegment}
end


