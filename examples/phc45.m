function p = phc45(N, L, Lframe, rhi, rlo, beta, Nz, betaz)
% 45 degree rotation of hexagonal phc, so that dipole x and y axes go along
% diagonals of grid. This is to address the issue that the TE Moperator has
% no mirror symmetry in x and y directions, but seems to have at least some
% degree of inversion symmetry. (It might have full inversion symmetry
% actually, but I'm not sure why the odd angular momentum TE disk modes,
% which seem to be fixed to the diagonals, are not exactly degenerate


% frame = 1 if you want to not stamp holes that get cut off at the edge of
% the rectangle.

    if ~exist('Nz')
        Nz = 1;
        betaz = 1;
    end

    h = L/(N-1);
    p = ones(N);
    
    xphole = [0, 1, 1, 1/2, 0];
    yphole = sqrt(3) * [1, 1, 0, 1/2, 0 ];
    
    for ix = 0:N-1,
        for iy = 0:N-1,
           
            x = (ix - (N-1)/2 ) * h;
            y = (iy - (N-1)/2 ) * h;
            
            xr = 1/sqrt(2) * ( x + y);
            yr = 1/sqrt(2) * ( -x + y);
            % rotated
            
            if abs(yr) > beta || abs(yr) > -sqrt(3)*abs(xr) + 2*beta,
                continue;
                % outside the hexagonal region
            end
            
            
            imx = floor( abs(xr) / 1 );
            imy = floor( abs(yr) / sqrt(3) );
            xp = mod(abs(xr), 1);
            yp = mod(abs(yr), sqrt(3) );
            % projected back into the first tile
            % floor and mod behave weirdly with negative numbers
            
            hole = 0; % records which hole, default 0 not in hole
            for ihole = 1:5,
                
                if( imx == 0 && imy == 0 && ihole == 5)
                    rcustom = rlo;
                else
                    rcustom = rhi;
                end
                
                if (xp - xphole(ihole))^2 + (yp - yphole(ihole))^2 < rcustom^2, 
                    hole = ihole;
                    break;
                end
            end
                        
            if hole == 0, 
                continue
            end
            
            xrhole = imx * 1 + xphole(hole); 
            yrhole = imy * sqrt(3) + yphole(hole);
            % note everything has already been absolute valued, so we are
            % only looking at upper right quadrant effectively
            
            xhole = 1/sqrt(2) * ( xrhole + -yrhole);
            yhole = 1/sqrt(2) * ( xrhole + yrhole);
            % rotate back  
            
            % make sure hole doesn't straddle edge of rectangular cell           
            if xhole + rhi < Lframe/2 && yhole + rhi < Lframe/2,
                p(iy+1, ix+1) = 0;
            end
        end 
    end

    % slab
    if Nz ~= 1
       p3d = zeros(Nz, N, N);        
       for iz=0:Nz-1
           z = (iz - (Nz-1)/2 ) * h;
           if abs(z) < betaz
              p3d(iz+1, :, :) = p;
           else
              p3d(iz+1, :, :) = zeros(size(p));
           end
       end
       p = p3d;        
    end


end
