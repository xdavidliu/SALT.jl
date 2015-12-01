N = 60; Nz = 1;
L = 6; Lframe = 100; 
format long;
h = L / (N-1)
Lz = (Nz-1)/(N-1) * L;
n_index = 11.56;
rhi = 0.3;
rlo = 0.2;
beta = 1.5; % for gain profile
betaz = 0.4; 
% if want no air, just periodic bc's in z-direction; i.e. pseudo-3d, then
% just set betaz to really high number
% if want slab, set betaz to less than Lz / 2

epsfile = 'eps3d.txt';
ffile = 'f3d.txt';

betaeps = 100; % basically infinite

p = phc45(N, L, Lframe, rhi, rlo, betaeps, Nz, betaz);
if Nz > 1, 
    pperm = permute(p, [2, 3, 1]);
    pslice = pperm(:, :, round(Nz/2));
else
    pslice = p;
end
figure
pcolor(pslice); axis equal; shading flat; axis off; 
title('eps')

eps = (n_index - 1)*p(:) + 1;
save(epsfile, 'eps', '-ascii', '-double');
Lframedummy = 1e5; 
% must be much greater than the large numbers used below for rhi and rlo when generating hexagonal region in fprof!

f = phc45(N, L, Lframedummy, 100, 100, beta, Nz, betaz);
fslab = phc45(N, L, Lframedummy, 0, 0, beta, Nz, betaz);
fprof = (1-f(:)) .* fslab(:);
% subtle: the inverting f with large rhi and rlo results in value of 1
% outside betaz. Pointwise multiply with an fslab with 0 outside betaz to
% obtain correct profile again

if Nz > 1,
    fperm = permute(f, [2, 3, 1]);
    fslice = 1-fperm(:, :, round(Nz/2));
else
    fslice = 1-f;
end
%figure;
%pcolor(fslice); axis equal; shading flat; axis off; 
%title('fprof');
% huge air holes to fill space inside beta

save(ffile, 'fprof', '-ascii', '-double');
