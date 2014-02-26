#include <slepc.h>


void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, char *epsfile, char *fproffile, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
);


int main(int argc, char** argv){

	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);


	int N[3] = {100, 1, 1}, Nc = 1,
	M[3] = {50, 1, 1},
	Npml[3] = {20, 0, 0},
	LowerPML = 0, bl[3] = {1, -1, 1},
	BCPeriod = -1, nev = 1,
	printnewton = 1
	;
	
	double h[3] = {0.01, 0.1, 0.2},
	wa = 15.0, y = 3.0, wreal = 15.7, wimag = -1.07,
	k[3] = {0., 0., 0.}, modenorm = 0.01,
	ftol = 1e-7, thresholdw_tol = 1e-7,
	dD = 0.05, Dmax = 0.0;
	;
	
	char epsfile[] = "eps1d.txt",
	fproffile[] = "fprof1d.txt",
	modeout[] = "pass1",
	Namesin[2][PETSC_MAX_PATH_LEN] = {"pass0", "pass1"},
	Namesout[2][PETSC_MAX_PATH_LEN] = {"last14", "last16"}
	; int Nm = 2;
	
	// hack because I don't know how const, char[][], char** and char *[] work
	char **namesin = (char**) malloc(2*sizeof(char*) ), 
	**namesout = (char**) malloc(2*sizeof(char*) );
		
	int i;
	for(i=0; i<2; i++){
		namesin[i] = (char*) malloc( PETSC_MAX_PATH_LEN *sizeof(char) );
		sprintf(namesin[i], "%s", Namesin[i]);
		namesout[i] = (char*) malloc( PETSC_MAX_PATH_LEN *sizeof(char) );
		sprintf(namesout[i], "%s", Namesout[i]);
	}

	Salt(N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y,
	BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout,
	dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm);

	SlepcFinalize();
}
