module SALT

export Passive, Creeper

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
		m_ = ccall( ("GetMode", saltlib), Mode_, (Ptr{Void}, Cint), ms, i-1)
		ma[i] = Mode( m_, geo);

		name = string("mode", string( ma[i].psi[end-1] )[1:4] ); # first few characters of real part of omega
		ccall( (:SetName, saltlib), Void, (Mode_, Ptr{Uint8}), m_, name);
		# just for the Newton solver print statements.
	end

	if Nadded==1
		return ma[1];
	else
		return ma;
	end
end

function Creeper(ms::Array{Mode, 1}, 
	geo::Geometry; Dmax::Cdouble=-5.0, ftol::Cdouble=1.0e-7, 
	printNewton::Bool=true, dD::Cdouble=0.0, steps::Int64=20)
# Dmax < 0 means stop at first threshold found.

	for i=1:length(ms)
		if ms[i].pump != ms[1].pump
			throw(ArgumentError("all modes must be at same pump strength"));
			return 0;
		end
	end

	if( Dmax <= ms[1].pump && Dmax >=0)
		throw(ArgumentError("Dmax must be higher than the Dinitial!"));
	end

	ccall( ("SetPump", saltlib), Void,
		(Geometry_, Cdouble), geo.geo, ms[1].pump); 

	if dD == 0.0
		dD = (Dmax - ms[1].pump) / steps;
		if(dD < 0) 
			dD = 0.05; # default value for step size
		end
	end

	Nm = size(ms, 1);
	msC = Array(Mode_, Nm );
	for i=1:Nm
		msC[i] = ccall( ("CopyMode", saltlib), Mode_, (Mode_,), ms[i].m ); 
	end

	Nlasing = ccall( ("Creeper", saltlib), Cint, 
		(Cdouble, Cdouble, Cdouble, Ptr{Mode_}, Cint, Cint, Geometry_),
		dD, Dmax, ftol, msC, int32(printNewton), int32(size(ms,1)), geo.geo  
	);
	print(Nlasing, " new modes lasing after Creeper\n");

	msout = Array(Mode, Nm);
	for i=1:Nm
		msout[i] = Mode( msC[i], geo );
	end
		
	msout;
end


end
