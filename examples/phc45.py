# module originally a matlab function

import numpy as np	

def phc45(N, L, Lframe, rhi, rlo, beta, Nz, betaz):


	h = L/(N-1.)
	p = np.ones([N,N])

	xphole = np.array( [0, 1, 1, 1./2, 0] )
	yphole = np.sqrt(3) * np.array( [1, 1, 0, 1./2, 0 ] )

	for ix in range(N):
		for iy in range(N):
	   
			x = (ix - (N-1.)/2 ) * h
			y = (iy - (N-1.)/2 ) * h

			xr = 1./np.sqrt(2) * ( x + y)
			yr = 1./np.sqrt(2) * ( -x + y)
			# rotated

			if np.abs(yr) > beta or np.abs(yr) > -np.sqrt(3)*np.abs(xr) + 2*beta: continue
			# outside the hexagonal region

		    
			imx = int( np.floor( np.abs(xr) / 1 ))
			imy = int( np.floor( np.abs(yr) / np.sqrt(3) ))
			xp = np.abs(xr) % 1
			yp = np.abs(yr) % np.sqrt(3) 
			# projected back into the first tile
			# floor and mod behave weirdly with negative numbers

			hole = -1 # records which hole, default -1 not in hole
			for ihole in range(5): # length of xphole and yphole

				if( imx == 0 and imy == 0 and ihole == 5-1):
					rcustom = rlo
				else:
					rcustom = rhi

				if (xp - xphole[ihole])**2 + (yp - yphole[ihole])**2 < rcustom**2: 
				    hole = ihole
				    break
				
			if hole == -1: continue

			xrhole = imx * 1 + xphole[hole] 
			yrhole = imy * np.sqrt(3) + yphole[hole]
			# note everything has already been absolute valued, so we are
			# only looking at upper right quadrant effectively

			xhole = 1/np.sqrt(2) * ( xrhole + -yrhole)
			yhole = 1/np.sqrt(2) * ( xrhole + yrhole)
			# rotate back  

			# make sure hole doesn't straddle edge of rectangular cell           
			if xhole + rhi < Lframe/2 and yhole + rhi < Lframe/2:
				p[iy, ix] = 0

	# slab
	if Nz != 1:
		p3d = np.zeros([Nz, N, N])        
		for iz in range(Nz):
			z =(iz - (Nz-1.)/2 ) * h
			if np.abs(z) < betaz:
				p3d[iz, :, :] = p
			else:
				p3d[iz, :, :] = np.zeros([len(p), len(p)])
		p = p3d

	return p        
