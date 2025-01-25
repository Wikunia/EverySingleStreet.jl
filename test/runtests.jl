using EverySingleStreet
using Test

using Dates
using Geodesy
using TimeZones
using LightXML

@testset "EverySingleStreet.jl" begin
    include("unit/download.jl")
    include("unit/gpx.jl")
    include("unit/candidates.jl")
    include("unit/map_path.jl")
    include("unit/shortest_path.jl")
    include("unit/districts.jl")
    include("unit/drawing.jl")
    include("unit/update_walked_parts.jl")
end
