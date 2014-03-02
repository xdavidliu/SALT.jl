#include "headers.h"

void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, double *eps, double *fprof, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
)
{
    Geometry geo;
    geo = CreateGeometry(N, M, h, Npml, Nc, LowerPML, eps, fprof, wa, y);    
    
    if(Dmax == 0.0){

		ModeArray ma = Passive(BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout, geo);
		int i;
		for(i = 0; i<ma->size; i++){
			Write(ma->L[i], geo);
			DestroyMode(ma->L[i]);
		}
		free(ma);
    }else
        Creeper(dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm, geo);    
    
    DestroyGeometry(geo);
}
