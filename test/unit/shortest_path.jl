@testset "get_shortest_path" begin 
    city_map = EverySingleStreet.parse_map(joinpath(@__DIR__, "..", "data", "CLZ.json"));
    ft = [(1517600713, 27031685), (5223629310, 30325936), (5223629310, 27031685), (5223629310, 1946216024)]
    for (from, to) in ft
        sp = EverySingleStreet.get_shortest_path(city_map, from, to)
        sp_old = EverySingleStreet.shortest_path(city_map.graph, from, to)
        @test sp == sp_old
    end
end