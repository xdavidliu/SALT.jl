using PyPlot

function quadrants( A, bx, by, centerstrip, hx, hy )

    A = [bx*fliplr(A) A[:, 1+centerstrip:end] ];  # note in Julia , also does vertical concatenation, so you need space here!
    A = [by*flipud(A); A[1+centerstrip:end, :] ];

    Nx = size(A, 2);
    Ny = size(A, 1);
    Lx = (Nx-1)*hx;
    Ly = (Ny-1)*hy;
    x = linspace(-Lx/2, Lx/2, Nx);
    y = linspace(-Ly/2, Ly/2, Ny)';

	X = ( ones(length(y))' .* x )'; # meshgrid using broadcasting
	Y = ( y .* ones(length(x)) )';
	
	return A, X, Y;
end

function plotTEslice(v, N, h, b, nz=1)

    Nx = N[1];
    Ny = N[2];
    Nz = N[3];
    
    hx = h[1];
    hy = h[2];
    bx = b[1];
    by = b[2];
            
    Nxyz = Nx * Ny * Nz;
    Exy = v[1:2*Nxyz] + im*v[3*Nxyz+1:5*Nxyz];
    misc = v[6*Nxyz+1:end];
 
    Ex = Exy[1:Nxyz];
    Ey = Exy[Nxyz+1:end];
    Ex = Ex[nz:Nz:Nxyz];
    Ey = Ey[nz:Nz:Nxyz]; # neat hack for slicing
    Exy = [Ex; Ey];
    v = [real(Exy); imag(Exy); misc];
    #================
    # from plotfield
    
    Ex = v[1:Nx*Ny] + im*v[2*Nx*Ny+1:3*Nx*Ny];
    Ex = reshape(Ex, Ny, Nx);
    Ey = v[Nx*Ny+1:2*Nx*Ny] + im*v[3*Nx*Ny+1:4*Nx*Ny];
    Ey = reshape(Ey, Ny, Nx);
    
    dyEx = diff(Ex,1)/hx;
    dxEy = diff(Ey,2)/hy; # julia's diff is missing the middle argument of Matlab's 3 arg version
    Hz = dxEy[1:end-1, :] - dyEx[:, 1:end-1];

    
    psi, X, Y = quadrants(Hz, -bx, -by, 0, hx, hy); 
        # magnetic field is pseudovector, so it works by the opposite BCs

    pcolor(X, Y, -real(psi), cmap="RdBu"); axis("equal"); axis("off");
end