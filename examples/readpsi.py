#!/usr/bin/python
filename = 'pass3d1_file.m'
Npml = [12, 12, 0]
import numpy as np

s = open(filename).read()
L = s.split("[\n")[1:]
# get rid of useless first three lines
K = L[0].split("\n];\n")
psi = [float(j) for j in K[0].split('\n')]
# get variables N, h, etc; only need stuff after k = ...
exec(L[1].split('\n];\n')[1])

Nxyz = N[0]*N[1]*N[2]
Nxyzc = Nxyz * Nc
Nr = 2
psicomp = np.array(psi[:Nxyzc]) + 1j*np.array(psi[Nxyzc:Nr*Nxyzc])
w = psi[-2]
c = 0
if psi[-1] < 0:
	w += psi[-1]*1j
else:
	c = psi[-1]

if Nc == 2 and N[2] == 1: # 2d TM
	Ex = np.reshape( psicomp[:Nxyz], [N[1], N[0]], 'F' )
	Ey = np.reshape( psicomp[Nxyz:2*Nxyz], [N[1], N[0]], 'F' )
	# ironically F means fortran order and is column major, 
	# which is also Xiangdong and Steven's convention
	# while C is C order and is row major
	dEydx = np.diff(Ey, axis=1)[:-1,:] / h[0]
	dExdy = np.diff(Ex, axis=0)[:,:-1] / h[1]
	curlE = dExdy - dEydx;
	curlEin = curlE[Npml[1]:-Npml[1], Npml[0]:-Npml[0]]

	from matplotlib import pyplot as plt
	plt.pcolor(np.real(curlEin), cmap='RdBu')
	plt.axis('equal')
	plt.savefig('pcolor.pdf')
