#include "headers.h"



void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, char *epsfile, char *fproffile, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
)
{

	Geometry Geo, *geo = &Geo;

	CreateGeometry(geo, N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);	
/*

	PetscPrintf(PETSC_COMM_WORLD, "DEBUG: Creating Geometry...\n");	
	PetscPrintf(PETSC_COMM_WORLD, "N = {%i, %i, %i}\n", N[0], N[1], N[2]);	
	PetscPrintf(PETSC_COMM_WORLD, "Npml = {%i, %i, %i}\n", Npml[0], Npml[1], Npml[2]);	
	PetscPrintf(PETSC_COMM_WORLD, "M = {%i, %i, %i}\n", M[0], M[1], M[2]);	
	PetscPrintf(PETSC_COMM_WORLD, "h = {%g, %g, %g}\n", h[0], h[1], h[2]);		
	PetscPrintf(PETSC_COMM_WORLD, "epsfile = %s\n", epsfile);
	PetscPrintf(PETSC_COMM_WORLD, "fproffile = %s\n", fproffile);		
	PetscPrintf(PETSC_COMM_WORLD, "wa = %g, y = %g\n", wa, y);





	PetscPrintf(PETSC_COMM_WORLD, "BCPeriod = %i, bl = {%i, %i, %i}, wreal = %g, wimag = %g\n",
		BCPeriod, bl[0], bl[1], bl[2], wreal, wimag);

	PetscPrintf(PETSC_COMM_WORLD, "k = {%g, %g, %g}, modenorm = %g\n", k[0], k[1], k[2], modenorm);
	 
	double epsnorm;
//	VecNorm(geo->veps, NORM_2, &epsnorm);

	VecSet(geo->vMscratch[0], 1.0);
	VecNorm(geo->vMscratch[0], NORM_2, &epsnorm);
	PetscPrintf(PETSC_COMM_WORLD, "|eps| = %g\n", epsnorm);

*/
	





	if(Dmax == 0.0)
		Passive(BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout, geo);		
	else
		Creeper(dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm, geo);		




	DestroyGeometry(geo);


}
