using LinearAlgebra

"""
    Line(n[, pt])

Store a line in N-dimensional space. The line passes through point `pt` with a normal direction `n`.
If `pt` is ommitted, it is taken to be the origin. If given, it must have the same length as `n`.
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
    HCylinder(radius, height)

Represent a cylinder of radius `radius` and `height`, oriented along the X axis.
"""
struct HCylinder <: Solid
    radius::Float64
    rad2::Float64
    height::Float64

    HCylinder(r, h) = new(r, r^2, h)
end

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
volume(c::HCylinder) = π*c.rad2*c.height
volume(b::Box) = prod(b.sides)

function path(box::Box, line::Line)
    return 0.0
end

"""
    path(obj::Solid, line::Line)

Return the path length of `line` through the object `obj`, or 0 if the line doesn't
intersect the object.
"""
function path(cyl::HCylinder, line::Line)
    # Step 0: If the line is parallel to the x-axis, special case.
    if line.n[2]==0 && line.n[3] == 0
        if line.pt[2]^2 + line.pt[3]^2 ≤ cyl.rad2
            return cyl.height
        end
        return 0.0
    end

    # Step 1: At which 0, 1, or 2 points does line intersect the infinitely long version of cyl?
    # Solving a quadratic equation
    a = line.n[2]^2+line.n[3]^2
    b = 2*(line.pt[2]*line.n[2] + line.pt[3]*line.n[3])
    c = line.pt[2]^2+line.pt[3]^2-cyl.rad2
    disc = b - 4a*c
    disc ≤ 0 && return 0 # negative or 0 discriminant means doesn't have finite-length intersection

    # Step 2: Solve quadratic. Find the points x1 and x2 where line enters/leaves the ∞ cylinder.
    avgt = -b/2a
    difft = sqrt(disc)/2a
    t1 = avgt-difft
    t2 = avgt+difft
    # Make sure that t1 has the lesser x component of the two. If not, swap.
    if t1[1] > t2[1]
        t1, t2 = t2, t1
    end

    # Find full points x1 and x2 (with x1[1] ≤ x2[1]) where line crosses the ∞ cylinder
    halfht = 0.5*cyl.height
    x1 = point(line, t1)
    x1[1] ≥ halfht && return 0.0
    x2 = point(line, t2)
    x2[1] ≤ -halfht && return 0.0

    # Step 3.1: if x1 is outside the finite cylinder, move it to point on the x=-h/2 boundary
    if x1[1] < -halfht
        t1 = -(line.pt[1]+halfht)/line.n[1]
    end
    # Step 3.2: if x2 is outside the finite cylinder, move it to point on the x=+h/2 boundary
    if x2[1] > +halfht
        t2 = -(line.pt[1]-halfht)/line.n[1]
    end
    abs(t2-t1)
end
