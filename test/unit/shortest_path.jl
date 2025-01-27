@testset "get_shortest_path" begin 
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    ft = [(1517600713, 27031685), (5223629310, 30325936), (5223629310, 27031685), (5223629310, 1946216024)]
    for (from, to) in ft
        sp = EverySingleStreet.get_shortest_path(city_map, from, to)
        sp_old = EverySingleStreet.shortest_path(city_map.graph, from, to)
        @test sp == sp_old
    end

    sp = EverySingleStreet.get_shortest_path(city_map, 30223832, 30223767; only_walkable_road=false)
    @test sp == [30223832, 30223767]

    sp = EverySingleStreet.get_shortest_path(city_map, 30223832, 30223767; only_walkable_road=true)
    @test sp == [ 30223832, 7147230001, 30337949, 30337959, 30337960, 7147229997, 30337961, 4734466864,30277127, 7147229999,2131474157, 30223767]
end