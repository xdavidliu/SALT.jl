immutable Mode_s; end
typealias Mode_ Ptr{Mode_s}
DestroyMode(m::Mode_) = ccall((:DestroyMode, saltlib), Void, (Mode_,), m )

function BCPeriod(N::Array{Int64,1}, k::Array{Cdouble,1})

	if( maximum( abs( k)) == 0.0)
		if( N[2] == 1 && N[3] == 1)
			return -1
		elseif( N[3] == 1)
			return 3
		else
			return 0
		end
	else
		throw(ArgumentError("automatic BCPeriod for nonzero k not implemented yet"));
		return 0;
	end
end

type Mode
	m::Mode_
	psi::PetscVec_
	N::Array{Int64,1}
	Nc::Int64
	h::Array{Cdouble,1}
	pump::Cdouble
	function Mode(m::Mode_, g::Geometry)
		md = new(
			m,
			ccall( (:GetVpsi, saltlib), PetscVec_, (Mode_,), m),
			GetN(g),
			GetNc(g),
			GetCellh(g),
			GetPump(g)
		);
		finalizer(md, DestroyMode)
		return md
	end
end

function Getbc(m::Mode)
	b = [0, 0, 0];
	for i=1:3
		b[i] = int64( ccall( (:Getbc, saltlib), Cint, (Mode_, Cint), m.m, i-1) );
	end
	b;
end

function show(io::IO, mode::Mode)
    print(io, "SALT Mode: ", mode.m, "\n")
	
	print(io, mode.Nc, " electric field components\n");
	print(io, mode.N[1], " x ", mode.N[2], " x ", mode.N[3], " pixels\n");
	h = mode.h;
	N = mode.N;
	Nc = mode.Nc;
	b = Getbc(mode);

	last_element = mode.psi[end];
	omega = mode.psi[end-1] + im * (last_element > 0? 0 : last_element); 
	magnitude = last_element > 0? last_element : 0;

	x = linspace(0, N[1]*h[1], N[1]);
	y = linspace(0, N[2]*h[2], N[2])';

	LowerPML = false; ## TODO
	figure();
	if( N[2]==1 && N[3] == 1 && Nc == 1)
		psi = mode.psi[1:N[1]] + im*mode.psi[N[1]+1:end-2];
		if !LowerPML
			x = [-flipud(x), x[2:end]];
			psi = [b[1]*flipud(psi), psi[2:end]]; 
		end	

		plot(x, real(psi), x, imag(psi));
		legend(["real", "imag"]);
		ylabel("\$\\Psi(x)\$");
		xlabel("\$x\$");
	elseif( N[2]>1 && Nc == 3)
		X = ( ones(length(y))' .* x )'; # meshgrid using broadcasting
		Y = ( y .* ones(length(x)) )';

		plotTEslice(mode.psi[1:end-2], N, h, b);
	end 

	title(
		string("Mode: \$\\omega\$ = ", string(real(omega))[1:min(7, end)], " + i(", 
		string(imag(omega))[1:min(7, end)], "), |\$\\Psi\$| = ", string(magnitude)[1:min(5, end)] )
	);

end