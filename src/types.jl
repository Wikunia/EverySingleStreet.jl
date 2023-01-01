
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
    nodes::Vector{Node}
    ways::Vector{Way}
end



