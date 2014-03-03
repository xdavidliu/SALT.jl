module SALT

export Salt, Passive, Creeper

#=========== libraries ============ #
const saltlib = Pkg.dir("SALT", "deps", "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], ENV["PETSC_ARCH"], "lib", "libslepc");

ccall((:SlepcInitialize, slepc), Cint,
	(Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
	C_NULL, C_NULL, C_NULL, C_NULL); # ignore arguments for now

###########################################################################
# Managing the Geometry type

immutable Geometry_s; end
typealias Geometry_ Ptr{Geometry_s}

immutable Mode_s; end
typealias Mode_ Ptr{Mode_s}


DestroyGeometry(geo::Geometry_) = ccall((:DestroyGeometry, saltlib), Void,
                                        (Geometry_,), geo)

DestroyMode(m::Mode_) = ccall((:DestroyMode, saltlib), Void, (Mode_,), m )

type Geometry
    geo::Geometry_
    function Geometry(geo::Geometry_)
        g = new(geo)
        finalizer(g, DestroyGeometry)
        return g
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

##################################################################

function Geometry(ε::Array{Cdouble}, h_, nPML_,
                  gain_prof::Array{Cdouble}, ω_gain::Real, γ_gain::Real; # keyword arguments next
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

##################################################################

function Passive(BCPeriod::Int64, bl::Array{Int64,1}, wreal::Cdouble, wimag::Cdouble, 
	geo::Geometry; nev::Integer=1, modenorm::Cdouble=0.1, k::Array{Cdouble,1} = [0.0, 0.0, 0.0])

	marray = ccall( ("Passive", saltlib), Ptr{Void}, (Cint, Ptr{Cint}, 
		Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, SALT.Geometry), 
		int32(BCPeriod), int32(bl), k, wreal, wimag, modenorm, nev, "", geo );

	N = ccall( ("GetArraySize", saltlib), Cint, (Ptr{Void},), marray );

	ma = Array(Mode, N);
	for i=1:N
		ma[i] = Mode( ccall( ("GetMode", saltlib), Mode_, (Ptr{Void}, Cint), marray, i-1) );
	end

	if nev==1
		return ma[1];
	else
		return ma;
	end
end

function Creeper(dD::Cdouble, Dinit::Cdouble, Dmax::Cdouble, ms::Array{Mode, 1}, 
	geo::Geometry; thresholdw_tol::Cdouble=1.0e-7, ftol::Cdouble=1.0e-7, 
	printNewton::Bool=true)

	
	ccall( ("SetPump", saltlib), Void,
		(Geometry_, Cdouble), geo.geo, Dinit); 


	msC = Array(Mode_, size(ms) );
	for i=1:size(ms,1)
		msC[i] = ms[i].m;
	end

	ccall( ("Creeper", saltlib), Void, 
		(Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Mode_}, Cint, Cint, Geometry_),
		dD, Dmax, thresholdw_tol, ftol, msC, int32(printNewton), int32(size(ms,1)), geo.geo  
	);

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




include("plotSalt.jl")

end
