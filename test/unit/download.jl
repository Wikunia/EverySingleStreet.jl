@testset "Download bounding box" begin 
    path = joinpath(@__DIR__, "..", "data", "test_download.json")
    EverySingleStreet.download(LLA(53.0, 9.9), LLA(53.01, 10.01) , path)
    city_map = EverySingleStreet.parse_map(path)
    @test !all(EverySingleStreet.iswalkable, city_map.ways)
    EverySingleStreet.filter_walkable_json!(path)
    city_map = EverySingleStreet.parse_map(path)
    @test all(EverySingleStreet.iswalkable, city_map.ways)
end

@testset "Download CLZ" begin 
    path = joinpath(@__DIR__, "..", "data", "CLZ.json")
    EverySingleStreet.download("Clausthal-Zellerfeld, Germany", path)
    city_map = EverySingleStreet.parse_map(path)
    @test !all(EverySingleStreet.iswalkable, city_map.ways)
    EverySingleStreet.filter_walkable_json!(path)
    city_map = EverySingleStreet.parse_map(path)
    @test all(EverySingleStreet.iswalkable, city_map.ways)
end