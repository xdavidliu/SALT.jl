iz = 7; % C 0-based index
Npmlx = 5;
Npmly = 5;

psicomp = psi(0+(1:prod(N)*3)) + 1i * psi((prod(N)*3)+(1:prod(N)*3));

E = {0,0,0}; for i=1:3, E{i} = psicomp((1:prod(N))+(i-1)*prod(N)); end
Eblock = {0,0,0}; for i=1:3, Eblock{i} = reshape(E{i}, N(3), N(2), N(1)); end
% column major order

Experm = permute(Eblock{1}, [2, 3, 1]);
Eyperm = permute(Eblock{2}, [2, 3, 1]);
% put z coordinate in very last so can broadcast x and y coordinates
% slightly confusing, because for column major, coordinates are z, y, x, so
% 1->z, 2->y, and 3->x

Exslice = Experm(:, :, iz+1);
Eyslice = Eyperm(:, :, iz+1);

% replace psi so plotrect_psiTE can work on it
psinew = [real(Exslice(:)); real(Eyslice(:)); imag(Exslice(:)); imag(Eyslice(:)); psi(end-1); psi(end)]; 

LowerPML = 1;
plotrect_psiTE(psinew, N(1), N(2), h(1), h(2), 1, 1, Npmlx, Npmly, LowerPML );