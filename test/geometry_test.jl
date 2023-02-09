using LinearAlgebra
using MuscRat

@testset "geometry constructors" begin
    l = Line([1,2,3], [4,5,6])
    @test norm(l.n) ≈ 1
    l4dim = Line([1,2,3,4], [4,5,6,7])
    @test norm(l4dim.n) ≈ 1
    @test_throws DimensionMismatch unequal_dim = Line([1,2,3,4], [1,2,3])

    hc = HCylinder(4, 5)
    @test hc.rad2 == 16.0
    @test volume(hc) ≈ 80π
    @test smallest_radius(hc) ≈ sqrt(4^2+2.5^2)

    box_vector = Box([3, 4, 5])
    box_tuple = Box(3, 4, 5.0)
    @test length(box_vector.sides) == 3
    @test length(box_tuple.sides) == 3
    @test volume(box_vector) ≈ 60
    @test volume(box_tuple) ≈ 60
    @test smallest_radius(box_vector) ≈ sqrt(1.5^2 + 2^2 + 2.5^2)
end

@testset "path lengths" begin
    @testset "HCylinder" begin
        R, H = 4, 10
        cyl = HCylinder(R, H)
        hline0 = Line([1, 0, 0])
        hlineout = Line([1, 0, 0], [0, R+1, 0])
        @test path(cyl, hline0) == H
        @test path(cyl, hlineout) == 0.0

        line_notouch = Line([0, 0, 1], [0, R+1, 0])
        line_tangent = Line([0, 0, 1], [0, R, 0])
        line_cross_left = Line([0, 0, 1], [-H, 0, 0])
        line_cross_right = Line([0, 0, 1], [+H, 0, 0])
        for line in (line_notouch, line_tangent, line_cross_left, line_cross_right)
            @test path(cyl,line) == 0.0
        end
        line_both_caps = Line([.9, 0, .1], [0, 0, 0])
        @test path(cyl, line_both_caps) ≈ H/cos(atan(1/9))
        line_no_caps = Line([1, 0, 1])
        @test path(cyl, line_no_caps) ≈ 2*sqrt(2)*R
        line_left_cap = Line([1, 0, 1], [-H/2, 0, 0])
        line_right_cap = Line([1, 0, 1], [H/2, 0, 0])
        @test path(cyl, line_left_cap) ≈ sqrt(2)*R
        @test path(cyl, line_right_cap) ≈ sqrt(2)*R
    end

    @testset "Box" begin
        b = Box([2,2,4])
        @test path(b, Line([1,1,0])) ≈ 2sqrt(2)
        @test path(b, Line([1,0,0])) ≈ 2
        @test path(b, Line([0,0,1])) ≈ 4
        path(b, Line([.5,1,0],[.5,0,0])) ≈ sqrt(5)
    end
end
