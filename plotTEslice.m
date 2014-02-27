function [X, Y, psi] = plotTEslice(v, N, h, b, nz)
% taken from slicefieldTE and plotfieldTE
 


    if(~exist('nz'))
        nz = 1;
    end

    Nx = N(1);
    Ny = N(2);
    Nz = N(3);
    
    hx = h(1);
    hy = h(2);
    bx = b(1);
    by = b(2);
    
        
    Nxyz = Nx * Ny * Nz;
 
    Exy = v(1:2*Nxyz) + 1i*v(3*Nxyz+1:5*Nxyz);
    misc = v(6*Nxyz+1:end);
 
    Ex = Exy(1:Nxyz);
    Ey = Exy(Nxyz+1:end);
    Ex = Ex(nz:Nz:Nxyz);
    Ey = Ey(nz:Nz:Nxyz); % neat hack for slicing
 
    Exy = [Ex; Ey];
 
    v = [real(Exy); imag(Exy); misc];
    
    %================
    % from plotfield
    
    Ex = v(1:Nx*Ny) + 1i*v(2*Nx*Ny+1:3*Nx*Ny);
    Ex = reshape(Ex, Ny, Nx);
    Ey = v(Nx*Ny+1:2*Nx*Ny) + 1i*v(3*Nx*Ny+1:4*Nx*Ny);
    Ey = reshape(Ey, Ny, Nx);
    
    dyEx = diff(Ex,1,1)/hx;
    dxEy = diff(Ey,1,2)/hy; % these are not square matrices
    Hz = dxEy(1:end-1, :) - dyEx(:, 1:end-1);
 
    
    [psi, X, Y] = quadrants(Hz, -bx, -by, 0, hx, hy); 
    %[psi, X, Y] = quadrants(Ex, -bx, by, 0, hx, hy); 
    %[psi, X, Y] = quadrants(Ey, bx, -by, 0, hx, hy); 
        % magnetic field is pseudovector, so it works by the opposite BCs
    bluered;
    pcolor(X, Y, -real(psi)); shading flat; axis equal; axis off; 
 
    X = X(:);
    Y = Y(:);
    psi = -real(psi(:));
    
end

function bluered(n)
  if nargin < 1
    n = 64;
  end
  r = linspace(0,1,n/2)';
  map = [r, r, ones(size(r))];
  map1 = flipud(fliplr(map(1:end-1,:)));
  map = [map; map1];
  colormap(map);
  c = max(abs(caxis));
  caxis([-c c]);
  
end