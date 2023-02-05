using EverySingleStreet
using Test

using Geodesy

@testset "EverySingleStreet.jl" begin
    include("unit/candidates.jl")
    include("unit/map_path.jl")
end
