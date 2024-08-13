
# double getHPcpp(int, int, int);
# double getrcpp(double, double);
# double getdcpp(double, double);
# double getSpecCpp(int, double, double, double, double, double);
# double getSpecAngFinalCpp(int, double, double, double, double, double, double);
# double get511fluxCpp(double, double, double);

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
    Parma_getHPcpp(year::Integer, month::Integer, day::Integer)

Return the W-index of solar activity at date `year`-`month`-`day`.
"""
function Parma_getHPcpp(year::Integer, month::Integer, day::Integer)
    depth = @ccall library._Z8getHPcppiii(year::Cint, month::Cint, day::Cint)::Float64
end
"""
    Parma_getdcpp(altitude::Real, latitude::Real=100)

Return the atmospheric depth at `alititude` (in km above sea level) and at 
`latitude` (in degrees), in g/cm^2. 

If |`latitude`| ≤ 90, use NRLMSISE-00 data.
Otherwise, use the U.S. Standard Atmosphere.
"""
function Parma_getdcpp(altitude::Real, latitude::Real=100)
    depth = @ccall library._Z7getdcppdd(altitude::Float64, latitude::Float64)::Float64
end