#!/usr/bin/python
N, Nz = 60, 25
L, Lframe = 6., 100. 
h = L / (N-1.)
print('h = %.15e' % h)
Lz = (Nz-1.)/(N-1) * L
n_index = 11.56
rhi = 0.3
rlo = 0.2
beta = 1.5 # for gain profile
betaz = 100 
# if want no air, just periodic bc's in z-direction i.e. pseudo-3d, then
# just set betaz to really high number
# if want slab, set betaz to less than Lz / 2

epsfilename = 'eps3d_py.txt'
ffilename = 'f3d_py.txt'

betaeps = 100 # basically infinite

import numpy as np
from phc45 import phc45
p = phc45(N, L, Lframe, rhi, rlo, betaeps, Nz, betaz)
if Nz > 1:
	pperm = p.transpose([1, 2, 0])
	pslice = pperm[:, :, Nz/2]
else: pslice = p

eps = (n_index - 1)*p + 1
epsfile = open(epsfilename, 'w')
epsC = eps.transpose([2, 1, 0]) # convert from fortran order to C order
epsC.tofile(epsfile, sep="\n", format="%.15e")
epsfile.close()

Lframedummy = 1e5 
# must be much greater than the large numbers used below for rhi and rlo when generating hexagonal region in fprof!

f = phc45(N, L, Lframedummy, 100, 100, beta, Nz, betaz)
fslab = phc45(N, L, Lframedummy, 0, 0, beta, Nz, betaz)
fprof = (1-f) * fslab
# subtle: the inverting f with large rhi and rlo results in value of 1
# outside betaz. Pointwise multiply with an fslab with 0 outside betaz to
# obtain correct profile again

if Nz > 1:
	fperm = f.transpose([1, 2, 0])
	fslice = 1-fperm[:, :, Nz/2]
else: fslice = 1-f

ffile = open(ffilename, 'w')
fprofC = fprof.transpose([2, 1, 0]) # convert from fortran order to C order
fprofC.tofile(ffile, sep="\n", format="%.15e")
ffile.close()


# must use ndarray.tofile(fid, sep="\n", format="%.15e")
# also this writes in C "major" order. Apparently matlab 3d arrays stored in fortran major order. Note that pF.transpose([2, 1, 0]) == pC. Hence, must do that transpose to the phc45 pblock array (which has the same exact elements in python and matlab, hence is in Fortran major order) before doing "tofile" on it
