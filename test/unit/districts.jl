@testset "get geojson file" begin
    districts = EverySingleStreet.get_districts(joinpath(@__DIR__, "..", "data", "hamburg_districts.geojson")); 
    @test length(districts) == 104
    district_names = [String(district.name) for district in districts]
    @test "Eimsbüttel" in district_names
    @test "Winterhude" in district_names
    @test "Bergedorf" in district_names
end

@testset "get district from node" begin
    nodes = [
        EverySingleStreet.Node(1, 53.57355, 9.95301) # in Eimsbüttel
        EverySingleStreet.Node(2, 53.5684, 9.9864) # in Rotherbaum
        EverySingleStreet.Node(3, 53.56155, 10.02219) # in Hohenfelde
    ]
    district_names = EverySingleStreet.map_nodes_to_district(nodes, joinpath(@__DIR__, "..", "data", "hamburg_districts.geojson"))
    @test district_names == [:Eimsbüttel, :Rotherbaum, :Hohenfelde]
end

