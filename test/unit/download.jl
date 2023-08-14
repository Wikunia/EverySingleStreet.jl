@testset "Download CLZ" begin 
    path = joinpath(@__DIR__, "..", "data", "CLZ.json")
    EverySingleStreet.download("Clausthal-Zellerfeld, Germany", path)
    city_map = EverySingleStreet.parse_map(path)
    @test !all(EverySingleStreet.iswalkable, city_map.ways)
    EverySingleStreet.filter_walkable_json!(path)
    city_map = EverySingleStreet.parse_map(path)
    @test all(EverySingleStreet.iswalkable, city_map.ways)
end