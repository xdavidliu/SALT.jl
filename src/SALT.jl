module SALT

export Salt, Passive, Creeper

#=========== libraries ============ #
const saltlib = Pkg.dir("SALT", "deps", "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], ENV["PETSC_ARCH"], "lib", "libslepc");

typealias PetscErrorCode Cint;
const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,slepc), Ptr{Void}));
ccall((:SlepcInitialize, slepc), PetscErrorCode,
					  (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              C_NULL, C_NULL, C_NULL, C_NULL);

###########################################################################
# Managing the Geometry type

immutable Geometry_s; end
typealias Geometry_ Ptr{Geometry_s}

immutable Mode_s; end
typealias Mode_ Ptr{Mode_s}

immutable ModeArray_s; end
typealias ModeArray_ Ptr{ModeArray_s}

DestroyGeometry(geo::Geometry_) = ccall((:DestroyGeometry, saltlib), Void,
                                        (Geometry_,), geo)

DestroyMode(m::Mode_) = ccall((:DestroyMode, saltlib), Void, (Mode_,), m )

DestroyModeArray(ma::ModeArray_) = ccall((:DestroyModeArray, saltlib), Void, (ModeArray_,), ma);


type Geometry
    geo::Geometry_
    function Geometry(geo::Geometry_)
        g = new(geo)
        finalizer(g, DestroyGeometry)
        return g
    end
end



type ModeArray
	ma::ModeArray_
	function ModeArray(ma::ModeArray_)
		mdar = new(ma)
		finalizer(mdar, DestroyModeArray)
		return mdar
	end
end


type Mode
	m::Mode_
	function Mode(m::Mode_)
		md = new(m)
		finalizer(md, DestroyMode)
		return md
	end


end


function Geometry(ε::Array{Cdouble}, h_, nPML_,
                  gain_prof::Array{Cdouble}, ω_gain::Real, γ_gain::Real;
                  nc::Integer=3, lowerPML::Bool=true)
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
    nPML = Cint[nPML_...]
    length(nPML) > 3 && throw(ArgumentError("nPML must have length <= 3"))
    while length(nPML) < 3
        push!(nPML, nPML[end])
    end


    Geometry(ccall(("CreateGeometry",saltlib), Geometry_,
                   (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
                    Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble),
                   N, N, h, nPML, nc, lowerPML, ε, gain_prof, ω_gain, γ_gain))
end



function ModeArray(BCPeriod::Int64, bl::Array{Int64,1}, k::Array{Cdouble,1},
    wreal::Cdouble, wimag::Cdouble, modenorm::Cdouble, geo::SALT.Geometry, n::Int64)


	ModeArray( ccall( ("Passive", saltlib), ModeArray_, (Cint, Ptr{Cint}, 
		Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, SALT.Geometry), 
		int32(BCPeriod), int32(bl), k, wreal, wimag, modenorm, 1, "", geo );
	)
end


function Mode(ma::ModeArray, n::Integer)
	Mode( ccall( ("GetMode", saltlib), Mode_, (ModeArray_, Cint), ma.ma, n) );
end



function GetPsi(m::SALT.Mode)
    
    N = ccall( (:PsiSize, saltlib), Cint, (Mode_,), m.m );
    v = zeros(N);
    ccall( (:CopyPsi, saltlib), Void, (Mode_, Ptr{Cdouble}), m.m, v);
    v;
end



import Base.show

function show(io::IO, g::Geometry)
    print(io, "SALT Geometry: ", g.geo)
end

###########################################################################

#=========== master Salt function from C code ============ #

Salt(N, nc, M, nPML, lowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm, nev) = 
ccall((:Salt, saltlib), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, 
Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble,
Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, 
Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Cint, 
Cint), int32(N), int32(M), h, int32(nPML), int32(nc), int32(lowerPML), 
eps, fprof, wa, y, int32(BCPeriod), int32(bl), k, wreal, wimag, modenorm, 
int32(nev), modeout, dD, Dmax, thresholdw_tol, ftol, namesin, namesout, 
int32(printnewton), int32(Nm));


#============ sub-function for Passive resonance eigensolver ========== #

Passive(N, nc, M, nPML, lowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout, nev) =
 
Salt(N, nc, M, nPML, lowerPML, bl, BCPeriod, 0, 
h, 0.0, 0.0, wreal, wimag, k, modenorm, 1.0e-7, 1.0e-7,
0.0, 0.0, eps, fprof, modeout, [""], [""], 0, nev);

#============ sub-function for Newton solver ======= #

Creeper(N, nc, M, nPML, lowerPML, printnewton, 
h, wa, y, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, namesin, namesout, Nm) =

Salt(N, nc, M, nPML, lowerPML, [0,0,0], 0, printnewton, 
h, wa, y, 0, 0, [0, 0, 0], 1.0, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, "", namesin, namesout, Nm, 0);

include("plotSalt.jl")

end
