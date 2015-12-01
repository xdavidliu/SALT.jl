function [X, Y, psi, Ex, Ey] = plotrect_psiTE(v, Nx, Ny, hx, hy, bx, by, Npmlx, Npmly, LowerPML )

	Ex = v(1:Nx*Ny) + 1i*v(2*Nx*Ny+1:3*Nx*Ny);
	Ex = reshape(Ex, Ny, Nx);
	Ey = v(Nx*Ny+1:2*Nx*Ny) + 1i*v(3*Nx*Ny+1:4*Nx*Ny);
	Ey = reshape(Ey, Ny, Nx);
	
	dyEx = diff(Ex,1,1)/hx;
	dxEy = diff(Ey,1,2)/hy; % these are not square matrices
	Hz = dxEy(1:end-1, :) - dyEx(:, 1:end-1);

	if(~exist('LowerPML', 'var')), LowerPML = 0; end
    
    if( LowerPML == 0),    
        [psi, X, Y] = quadrants(Hz, -bx, -by, 0, hx, hy); 
	%[psi, X, Y] = quadrants(Ex, -bx, by, 0, hx, hy); 
	%[psi, X, Y] = quadrants(Ey, bx, -by, 0, hx, hy); 
		% magnetic field is pseudovector, so it works by the opposite BCs
        
        ind_nopmlx = Npmlx+1:(2*Nx-2)-Npmlx; % need to subtract 2 because of funny business with Hz missing a row/column
        ind_nopmly = Npmly+1:(2*Ny-2)-Npmly;
        
    else
        psi = Hz;
        [X, Y] = meshgrid(0:hx:(Nx-2)*hx, 0:hy:(Ny-2)*hy); % should be -1, but Hz has one less row and column

        ind_nopmlx = Npmlx+1:Nx-1-Npmlx; % need to subtract 1 because of funny business with Hz missing a row/column
        ind_nopmly = Npmly+1:Ny-1-Npmly;
	
    end
    X = X(ind_nopmly, ind_nopmlx);
    Y = Y(ind_nopmly, ind_nopmlx);
    psi = psi(ind_nopmly, ind_nopmlx);
    
	pcolor(X, Y, real(psi)); shading flat; axis equal; axis off; bluered;	   

end

