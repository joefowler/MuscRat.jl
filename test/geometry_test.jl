using LinearAlgebra
using MuscRat

@testset "geometry constructors" begin
    l = Line([1,2,3], [4,5,6])
    @test norm(l.n) ≈ 1
    l4 = Line([1,2,3,4], [4,5,6,7])
    @test norm(l4.n) ≈ 1
    @test_throws DimensionMismatch unequal_dim = Line([1,2,3,4], [1,2,3])

    hc = HCylinder(4, 5)
    @test hc.rad2 == 16.0
    @test volume(hc) ≈ 80π

    box_vector = Box([3, 4, .5])
    box_tuple = Box(3, 4, .5)
    @test length(box_vector.sides) == 3
    @test length(box_tuple.sides) == 3
    @test volume(box_vector) ≈ 6
    @test volume(box_tuple) ≈ 6
end
