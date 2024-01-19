using EverySingleStreet
using Test

using Dates
using Geodesy
using TimeZones

@testset "EverySingleStreet.jl" begin
    include("unit/utils.jl")
    include("unit/download.jl")
    include("unit/gpx.jl")
    include("unit/candidates.jl")
    include("unit/map_path.jl")
    include("unit/shortest_path.jl")
    include("unit/districts.jl")
end
