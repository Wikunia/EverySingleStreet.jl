using EverySingleStreet
using Test

using Dates
using Geodesy
using TimeZones

@testset "EverySingleStreet.jl" begin
    include("unit/gpx.jl")
    include("unit/candidates.jl")
    include("unit/map_path.jl")
end
