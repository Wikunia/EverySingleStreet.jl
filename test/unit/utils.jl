@testset "config" begin 
    @test parse(Int, EverySingleStreet.CONFIG["GPS_STD_DEV"]) == 10
end