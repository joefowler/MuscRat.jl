using LinearAlgebra

"""
    Line(n[, pt])

Store a line in N-dimensional space. The line passes through point `pt` with a normal direction `n`.
If `pt` is omitted, it is taken to be the origin. If given, it must have the same length as `n`.
"""
struct Line
    n::Vector{Float64}
    pt::Vector{Float64}

    function Line(n::AbstractVector{T}, pt::AbstractVector{U}) where {T,U <: Real}
        if length(n) != length(pt)
            throw(DimensionMismatch("Direction vector and point must have equal length"))
        end
        new(n/norm(n), pt)
    end
end
Line(n::AbstractVector{T}) where {T<:Real} = Line(n, zero(n))
point(line::Line, t::Real) = line.n*t+line.pt

abstract type Solid end

"""
    Sphere(radius)

Represent a sphere of radius `radius`.
"""
struct Sphere <: Solid
    radius::Float64
    rad2::Float64

    Sphere(r) = new(r, r^2)
end

"""
    Cylinder(radius, height, axis)

Represent a cylinder of radius `radius` and `height`, oriented along the `axis` axis (must be 1, 2, or 3).
"""
struct Cylinder <: Solid
    radius::Float64
    rad2::Float64
    height::Float64
    axis::Int

    function Cylinder(r::Real, h::Real, axis::Integer)
        if !(axis in (1,2,3))
            error("Cylinder axis must be one of {1,2,3}")
        end
        new(r, r^2, h, axis)
    end
end
HCylinder(r::Real, h::Real) = Cylinder(r, h, 1)
VCylinder(r::Real, h::Real) = Cylinder(r, h, 3)

"""
    Box(sides)
    Box(a, b, c)

Represent an N-dimensional rectangular prism centered on the origin.
"""
struct Box <: Solid
    sides::Vector{Float64}
end
Box(args...) = Box([args...])

"""
    volume(obj::Solid)

Returns the volume of a solid object
"""
volume(s::Sphere) = (4/3)*π*(s.radius^3)
volume(c::Cylinder) = π*c.rad2*c.height
volume(b::Box) = prod(b.sides)

"""
    smallest_radius(obj::Solid)

Returns the radius of the smallest cylinder that, when centered on the origin at an arbitrary orientation,
will contain every point of the object `obj`.
"""
smallest_radius(s::Sphere) = s.radius
smallest_radius(c::Cylinder) = sqrt(c.radius^2+(c.height*0.5)^2)
smallest_radius(b::Box) = norm(b.sides)/2

"""
    tube_area(obj::Solid)

Returns the cross-sectional area of the smallest-radius tube that, if centered on `obj`, will
always enclose the object, regardless of the tube orientation.
"""
tube_area(obj::Solid) = π*smallest_radius(obj)^2

function path(box::Box, line::Line)
    crossings = []
    for i=1:3
        line.n[i]==0 && continue
        for s in [-0.5, +0.5]*box.sides[i]
            t = (s-line.pt[i])/line.n[i]
            x = point(line, t)
            inrange = true
            for j=1:3
                j==i && continue
                if abs(x[j]) > box.sides[j]*0.5
                    inrange = false
                    break
                end
            end
            if inrange
                push!(crossings, x)
            end
        end
    end
    length(crossings) == 0 && return 0.0
    # There should be exactly 2 good points.
    # The same 2 points can appear 2x each, though, if they are at the edges or corners of the box.
    # Even in that "corner case", the first 2 are guaranteed to be opposite one another.
    norm(crossings[1]-crossings[2])
end

"""
    path(obj::Solid, line::Line)

Return the path length of `line` through the object `obj`, or 0 if the line doesn't
intersect the object.
"""
function path(cyl::Cylinder, line::Line)
    # Step 1: At which 0, 1, or 2 points does line intersect the infinitely long version of cyl?
    # This means solving a quadratic equation (in general, though not quadratic if a=0).
    # r2 is the square distance of the line.pt from the cylinder axis
    a = b = r2 = 0.0
    for i=1:3
        if i != cyl.axis
            a += line.n[i]^2
            b += 2*line.pt[i]*line.n[i]
            r2 += line.pt[i]^2
        end
    end
    c = r2-cyl.rad2
    if a == 0
        if c ≤ 0
            return cyl.height
        end
        return 0.0
    end
    disc = b - 4a*c
    disc ≤ 0 && return 0 # negative or 0 discriminant means there's no intersection, or a length-0 tangent point only

    # Step 2: Solve quadratic. Find the points x1 and x2 where line enters/leaves the ∞ cylinder.
    avgt = -b/2a
    difft = sqrt(disc)/2a
    t1 = avgt-difft
    t2 = avgt+difft
    x1 = point(line, t1)
    x2 = point(line, t2)
    # Make sure that x1 is the lesser of the two along the cylinder axis. If not, swap.
    if x1[cyl.axis] > x2[cyl.axis]
        t1, t2 = t2, t1
        x1, x2 = x2, x1
    end

    # Find full points x1 and x2 (with x1[cyl.axis] ≤ x2[cyl.axis]) where line crosses the ∞ cylinder
    halfht = 0.5*cyl.height
    x1[cyl.axis] ≥ halfht && return 0.0
    x2[cyl.axis] ≤ -halfht && return 0.0

    # Step 3.1: if x1 is outside the finite cylinder, move it to point on the x=-h/2 boundary
    if x1[cyl.axis] < -halfht
        t1 = -(line.pt[cyl.axis]+halfht)/line.n[cyl.axis]
    end
    # Step 3.2: if x2 is outside the finite cylinder, move it to point on the x=+h/2 boundary
    if x2[cyl.axis] > +halfht
        t2 = -(line.pt[cyl.axis]-halfht)/line.n[cyl.axis]
    end
    abs(t2-t1)
end

function path(s::Sphere, line::Line)
    # Find the points (if any) where t allows ||pt+n⋅t|| = r^2
    # a = 1.0
    b = 2dot(line.n, line.pt)
    c = dot(line.pt, line.pt) - s.rad2
    disc = b - 4c
    disc ≤ 0 && return 0 # negative or 0 discriminant means doesn't have finite-length intersection

    # avgt = -b/2
    # difft = sqrt(disc)/2
    # t1 = avgt-difft
    # t2 = avgt+difft

    # The distance between these points is the abs difference of the t values
    sqrt(disc)
end