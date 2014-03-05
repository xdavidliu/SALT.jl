immutable Geometry_s; end
typealias Geometry_ Ptr{Geometry_s}
DestroyGeometry(geo::Geometry_) = ccall((:DestroyGeometry, saltlib), Void,
                                        (Geometry_,), geo)

type Geometry
    geo::Geometry_
	eps::PetscVec_
	fprof::PetscVec_
    function Geometry(geo::Geometry_)
        g = new(
			geo,
			ccall( (:GetVeps, saltlib), PetscVec_, (Geometry_,), geo),
			ccall( (:GetVfprof, saltlib), PetscVec_, (Geometry_,), geo)		
		)
        finalizer(g, DestroyGeometry)
        return g
    end
end

function Geometry(ε::Array{Cdouble}, grid_spacing, nPML_,
                  gain_prof::Array{Cdouble}, ω_gain::Real, γ_gain::Real;
                  n_vectorcomp::Integer=3, lowerPML::Bool=true)
    ndims(ε) > 3 && throw(ArgumentError("ε array must be <= 3-dimensional"))
    size(ε) != size(gain_prof) && throw(ArgumentError("gain profile must be same size as ε array"))
    N = Cint[size(ε)...]
	N = flipud(N); #ε column major, want N = Nx, Ny, Nz 

    while length(N) < 3
        push!(N, 1)
    end

    h = Cdouble[grid_spacing...]
    length(h) > 3 && throw(ArgumentError("h must have length <= 3"))
    while length(h) < 3
        push!(h, h[end])
    end
    nPML = Cint[nPML_...]
    length(nPML) > 3 && throw(ArgumentError("nPML must have length <= 3"))
    while length(nPML) < 3
        push!(nPML, nPML[end])
    end

    Geometry( ccall( ("CreateGeometry",saltlib), Geometry_,
            (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
            Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble),
            N, h, nPML, n_vectorcomp, lowerPML, ε, gain_prof, ω_gain, γ_gain))
end

# TODO: clear up names of Julia type and C objects, i.e. no geo.geo
function GetN(geo::Geometry)
	N = [0, 0, 0];
	for i=1:3
		N[i] = int64( ccall( (:GetN, saltlib), Cint, (Geometry_, Cint), geo.geo, i-1) );
	end
	N;
end

function GetNc(geo::Geometry)
	ccall( (:GetNc, saltlib), Cint, (Geometry_,), geo.geo);
end


function GetNpml(geo::Geometry)
	N = [0, 0, 0];
	for i=1:3
		N[i] = int64( ccall( (:GetNpml, saltlib), Cint, (Geometry_, Cint), geo.geo, i-1) );
	end
	N;
end

function GetCellh(geo::Geometry)
	h = [0.0, 0.0, 0.0];
	for i=1:3
		h[i] = ccall( (:GetCellh, saltlib), Cdouble, (Geometry_, Cint), geo.geo, i-1);
	end
	h;
end

function GetPump(geo::Geometry)
	ccall( (:GetD, saltlib), Cdouble, (Geometry_,), geo.geo)
end

function show(io::IO, g::Geometry)
    print(io, "SALT Geometry: ", g.geo, "\n")
	
	Nc = GetNc(g);

	N = GetN(g);
	print(io, N[1], " x ", N[2], " x ", N[3], " pixels in cell\n");
	Npml = GetNpml(g);
	print(io, Npml[1], " x ", Npml[2], " x ", Npml[3], " pixel PML thickness\n");
	h = GetCellh(g);
	print(io, h[1], " x ", h[2], " x ", h[3], " cell\n");
	x = linspace(0, N[1]*h[1], N[1]);	
	y = linspace(0, N[2]*h[2], N[2])';
	eps = g.eps[1:N[1]*N[2]*N[3]]; # real part only
	fprof = g.fprof[1:N[1]*N[2]*N[3]];

	subplot(221)

	# cannot get L"$\epsilon" to work here
	if( N[2]==1 && N[3] == 1)
		plot(x, eps);
		xlabel("x");
		ylabel("\$\\epsilon\$");
	elseif( N[2]>1 && N[3]>1 )
		X = ( ones(length(y))' .* x )'; # meshgrid using broadcasting
		Y = ( y .* ones(length(x)) )';
		eps = reshape(eps, N[3], N[2], N[1]);

		pcolor(X, Y, squeeze( eps[1, :, :],1 ), cmap="Blues" );
		axis("equal"); axis("off");
		colorbar();
	end
	title("Dielectric of passive medium");

	subplot(222)
	if( N[2]==1 && N[3] == 1)
		plot(linspace(0, N[1]*h[1], N[1]), fprof);
		xlabel("x");
		ylabel("\$f(x)\$");	
	elseif( N[2]>1 && N[3]>1 )
		fprof = reshape(fprof, N[3], N[2], N[1]);
		X = ( ones(length(y))' .* x )'; # meshgrid using broadcasting
		Y = ( y .* ones(length(x)) )';
		pcolor(X, Y, squeeze( fprof[1, :, :],1 ), cmap="Greens" );
		axis("equal"); axis("off");		
	end
	title("gain profile");

end