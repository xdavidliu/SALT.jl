module SALT

export Salt, Passive, Creeper

#=========== libraries ============ #
const saltlib = Pkg.dir("SALT", "deps", "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], get(ENV, "PETSC_ARCH", ""), "lib", "libslepc");

typealias PetscErrorCode Cint;
const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,slepc), Ptr{Void}));
ccall((:SlepcInitialize, slepc), PetscErrorCode,
              (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              C_NULL, C_NULL, C_NULL, C_NULL);

###########################################################################
# Managing the Geometry type

immutable Geometry_s; end
typealias Geometry_ Ptr{Geometry_s}

DestroyGeometry(geo::Geometry_) = ccall((:DestroyGeometry, saltlib), Void,
                                        (Geometry_,), geo)

type Geometry
    geo::Geometry_
    function Geometry(geo::Geometry_)
        g = new(geo)
        finalizer(g, DestroyGeometry)
        return g
    end
end

function Geometry(ε::Array{Cdouble}, h_, Npml_,
                  gain_prof::Array{Cdouble}, ω_gain::Real, γ_gain::Real;
                  Nc::Integer=3, LowerPML::Bool=true)
    ndims(ε) > 3 && throw(ArgumentError("ε array must be <= 3-dimensional"))
    size(ε) != size(gain_prof) && throw(ArgumentError("gain profile must be same size as ε array"))
    N = Cint[size(ε)...]
    while length(N) < 3
        push!(N, 1)
    end
    h = Cdouble[h_...]
    length(h) > 3 && throw(ArgumentError("h must have length <= 3"))
    while length(h) < 3
        push!(h, h[end])
    end
    Npml = Cint[Npml_...]
    length(Npml) > 3 && throw(ArgumentError("Npml must have length <= 3"))
    while length(Npml) < 3
        push!(Npml, Npml[end])
    end

    Geometry(ccall(("CreateGeometry",saltlib), Geometry_,
                   (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
                    Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble),
                   N, N, h, Npml, Nc, LowerPML, ε, gain_prof, ω_gain, γ_gain))
end

import Base.show

function show(io::IO, g::Geometry)
    print(io, "SALT Geometry: ", g.geo)
end

###########################################################################

#=========== master Salt function from C code ============ #

Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm) = 
ccall((:Salt, saltlib), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, 
Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble,
Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, 
Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Cint, 
Cint), int32(N), int32(M), h, int32(Npml), int32(Nc), int32(LowerPML), 
eps, fprof, wa, y, int32(BCPeriod), int32(bl), k, wreal, wimag, modenorm, 
int32(nev), modeout, dD, Dmax, thresholdw_tol, ftol, namesin, namesout, 
int32(printnewton), int32(Nm));


#============ sub-function for Passive resonance eigensolver ========== #

Passive(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout) =
 
Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 0, 
h, 0.0, 0.0, wreal, wimag, k, modenorm, 1.0e-7, 1.0e-7,
0.0, 0.0, eps, fprof, modeout, [""], [""], 0);

#============ sub-function for Newton solver ======= #

Creeper(N, Nc, M, Npml, LowerPML, printnewton, 
h, wa, y, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, namesin, namesout, Nm) =

Salt(N, Nc, M, Npml, LowerPML, [0,0,0], 0, printnewton, 
h, wa, y, 0, 0, [0, 0, 0], 1.0, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, "", namesin, namesout, Nm);

include("plotSalt.jl")

end
