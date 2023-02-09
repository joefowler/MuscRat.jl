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

function path(cyl::HCylinder, line::Line)
    return 0.0
end