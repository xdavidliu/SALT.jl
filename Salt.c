#include "headers.h"



void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, char *epsfile, char *fproffile, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
)
{


	Vec veps, vfprof;
	CreateVec(M[0]*M[1]*M[2], &veps);
	VecDuplicate(veps, &vfprof);

	FILE *fp;
	
	fp = fopen(epsfile, "r");
	if(fp==NULL){
		char message[PETSC_MAX_PATH_LEN];
		sprintf(message, "failed to read %s", epsfile);
		MyError(message);
	}
	ReadVectorC(fp, M[0]*M[1]*M[2], veps);
	fclose(fp);

	fp = fopen(fproffile, "r");	
	if(fp==NULL){
		char message[PETSC_MAX_PATH_LEN];
		sprintf(message, "failed to read %s", fproffile);
		MyError(message);
	}	
	ReadVectorC(fp, M[0]*M[1]*M[2], vfprof);
	fclose(fp);	





	double *eps, *fprof;
	VecGetArray(veps, &eps);
	VecGetArray(vfprof, &fprof);

	Geometry Geo, *geo = &Geo;
	CreateGeometry(geo, N, M, h, Npml, Nc, LowerPML, eps, fprof, wa, y);	

	VecRestoreArray(veps, &eps);
	VecRestoreArray(vfprof, &fprof);

	VecDestroy(&veps);
	VecDestroy(&vfprof);
	

	if(Dmax == 0.0)
		Passive(BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout, geo);		
	else
		Creeper(dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm, geo);		




	DestroyGeometry(geo);


}
