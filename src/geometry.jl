using LinearAlgebra

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

struct HCylinder
    radius::Float64
    rad2::Float64
    height::Float64

    HCylinder(r, h) = new(r, r^2, h)
end

struct Box
    sides::Vector{Float64}
end
Box(args...) = Box([args...])

function path(box::Box, line::Line)
    return 0.0
end

function path(cyl::HCylinder, line::Line)
    return 0.0
end