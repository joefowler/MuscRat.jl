using LinearAlgebra
using MuscRat

@testset "geometry constructors" begin
    l = Line([1,2,3], [4,5,6])
    @test norm(l.n) ≈ 1
    l4dim = Line([1,2,3,4], [4,5,6,7])
    @test norm(l4dim.n) ≈ 1
    @test_throws DimensionMismatch unequal_dim = Line([1,2,3,4], [1,2,3])

    hc = HCylinder(4, 5)
    vc = VCylinder(4, 5)
    c2 = Cylinder(4, 5, 2)
    for cyl in (hc, vc, c2)
        @test cyl.rad2 == 16.0
        @test volume(cyl) ≈ 80π
        @test smallest_radius(cyl) ≈ sqrt(4^2+2.5^2)
    end
    @test_throws MethodError Cylinder(4, 5)
    @test_throws ErrorException Cylinder(4, 5, 93)

    s = Sphere(4)
    @test s.rad2 == 16.0
    @test volume(s) ≈ 256π/3
    @test smallest_radius(s) ≈ 4.0

    box_vector = Box([3, 4, 5])
    box_tuple = Box(3, 4, 5.0)
    for box in (box_vector, box_tuple)
        @test length(box.sides) == 3
        @test volume(box) ≈ 60
        @test smallest_radius(box) ≈ sqrt(1.5^2 + 2^2 + 2.5^2)
    end
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
    @testset "VCylinder" begin
        R, H = 4, 10
        cyl = VCylinder(R, H)
        hline0 = Line([0, 0, 1])
        hlineout = Line([0, 0, 1], [0, R+1, 0])
        @test path(cyl, hline0) == H
        @test path(cyl, hlineout) == 0.0
        line_notouch = Line([1, 0, 0], [0, R+1, 0])
        line_tangent = Line([1, 0, 0], [0, R, 0])
        line_cross_below= Line([1, 0, 0], [0, 0, -H])
        line_cross_above = Line([1, 0, 0], [0, 0, +H])
        for line in (line_notouch, line_tangent, line_cross_below, line_cross_above)
            @test path(cyl,line) == 0.0
        end
        line_both_caps = Line([.3, 0, .8], [0, 0, 0])
        @test path(cyl, line_both_caps) ≈ H/cos(atan(3/8))
    end

    @testset "Box" begin
        b = Box([2,2,4])
        @test path(b, Line([1,1,0])) ≈ 2sqrt(2)
        @test path(b, Line([1,0,0])) ≈ 2
        @test path(b, Line([0,0,1])) ≈ 4
        @test path(b, Line([.5,1,0],[.5,0,0])) ≈ sqrt(5)
    end

    @testset "Sphere" begin
        s = Sphere(4)
        @test path(s, Line(randn(3))) ≈ 8
        for offset in 1:6
            p = path(s, Line([1,0,0], [0,offset,0]))
            if offset ≥ s.radius
                @test p == 0.0
            else
                @test p ≈ 2sqrt(s.rad2-offset^2)
            end
        end
    end
end
