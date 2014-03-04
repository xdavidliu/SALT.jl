module SALT

export Salt, Passive, Creeper

#=========== libraries ============ #
const saltlib = Pkg.dir("SALT", "deps", "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], ENV["PETSC_ARCH"], "lib", "libslepc");

ccall((:SlepcInitialize, slepc), Cint,
	(Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
	C_NULL, C_NULL, C_NULL, C_NULL); # ignore arguments for now

import Base.show

include("PetscVec.jl");
include("Geometry.jl");
include("Mode.jl");
include("plotSalt.jl");

##################################################################

function Passive(boundary_condition, ωguess, geo::Geometry; 
	nev::Integer=1, modenorm::Cdouble=0.1, k::Array{Cdouble,1} = [0.0, 0.0, 0.0])

    bl = Cint[boundary_condition...]
    while length(bl) < 3
        push!(bl, bl[end])
    end

	N = GetN(geo);
	Nadded = [0];
	ms = ccall( ("Passive", saltlib), Ptr{Void}, (Ptr{Cint}, Cint, Ptr{Cint}, 
		Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Geometry), 
		Nadded, int32( BCPeriod(N, k) ), int32(bl), k, real(ωguess), imag(ωguess), modenorm, nev, geo );
	Nadded = Nadded[1];

	ma = Array(Mode, Nadded);
	for i=1:Nadded
		ma[i] = Mode( 
			ccall( ("GetMode", saltlib), Mode_, (Ptr{Void}, Cint), ms, i-1), 
			geo
		);
	end

	if Nadded==1
		return ma[1];
	else
		return ma;
	end
end

function Creeper(dD::Cdouble, Dinit::Cdouble, Dmax::Cdouble, ms::Array{Mode, 1}, 
	geo::Geometry; thresholdw_tol::Cdouble=1.0e-7, ftol::Cdouble=1.0e-7, 
	printNewton::Bool=true)

	Nm = size(ms, 1);
	
	ccall( ("SetPump", saltlib), Void,
		(Geometry_, Cdouble), geo.geo, Dinit); 


	msC = Array(Mode_, Nm );
	for i=1:Nm
		msC[i] = ccall( ("CopyMode", saltlib), Mode_, (Mode_,), ms[i].m ); 
	end

	ccall( ("Creeper", saltlib), Void, 
		(Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Mode_}, Cint, Cint, Geometry_),
		dD, Dmax, thresholdw_tol, ftol, msC, int32(printNewton), int32(size(ms,1)), geo.geo  
	);

	msout = Array(Mode, Nm);
	for i=1:Nm
		msout[i] = Mode( msC[i], geo );
	end
		
	msout;
end


end
