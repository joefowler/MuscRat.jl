
# double getSpecAngFinalCpp(int, double, double, double, double, double, double);
# double get511fluxCpp(double, double, double);
# __Z18getSpecAngFinalCppidddddd
# __Z13get511fluxCppddd

@enum ParmaParticle begin
    ppNeutron = 0
    ppProton = 1
    ppHelium = 2
    ppµplus = 29
    ppµminus = 30
    ppElectron = 31
    ppPositron = 32
    ppGamma = 33
end

library = "parma/parma"

"""
    Parma_getHP(year::Integer, month::Integer, day::Integer)

Return the W-index of solar activity at date `year`-`month`-`day`.
"""
function Parma_getHP(year::Integer, month::Integer, day::Integer)
    Windex = @ccall library._Z8getHPcppiii(year::Cint, month::Cint, day::Cint)::Cdouble
end


"""
    Parma_getr(latitude::Real, longitude::Real=100)

Return the vertical cut-off rigidity (GV) at `latitude` and `longitude` (both in degrees).
"""
function Parma_getr(latitude::Real, longitude::Real=100)
    rigidity = @ccall library._Z7getrcppdd(latitude::Cdouble, longitude::Cdouble)::Cdouble
end

"""
    Parma_getd(altitude::Real, latitude::Real=100)

Return the atmospheric depth at `alititude` (in km above sea level) and at 
`latitude` (in degrees), in g/cm^2. 

If |`latitude`| ≤ 90, use NRLMSISE-00 data.
Otherwise, use the U.S. Standard Atmosphere.
"""
function Parma_getd(altitude::Real, latitude::Real=100)
    depth = @ccall library._Z7getdcppdd(altitude::Cdouble, latitude::Cdouble)::Cdouble
end

