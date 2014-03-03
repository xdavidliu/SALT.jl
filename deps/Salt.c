#include "headers.h"

void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, double *eps, double *fprof, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
)
{
    Geometry geo;
    geo = CreateGeometry(N, M, h, Npml, Nc, LowerPML, eps, fprof, wa, y);    
    
    if(Dmax == 0.0){

		ModeArray ma = Passive(BCPeriod, bl, k, wreal, wimag, modenorm, nev, geo);
		int i;
		for(i = 0; i<ma->size; i++){

			if(ma->size == 1)
				sprintf(ma->L[i]->name, "%s", modeout);
			else 
				sprintf(ma->L[i]->name, "%s%i", modeout, i);


			Write(ma->L[i], geo);
			DestroyMode(ma->L[i]);
		}
		free(ma);

		double *vals;
		VecGetArray(geo->vscratch[0], &vals);

		CopyPsi(ma->L[0], vals);
		VecRestoreArray(geo->vscratch[0], &vals);

		Output(geo->vscratch[0], "vals", "");

    }else{

	  	ModeArray ma = CreateModeArray();
		ReadModes(ma, geo, namesin, namesout, Nm);

		// hack, so that SALT.jl can access Creeper using just a pointer to modes.
        Creeper(dD, Dmax, thresholdw_tol, ftol, ma->L, printnewton, Nm, geo);    

		int ih;

		for(ih=0; ih<ma->size; ih++){
			Write(ma->L[ih], geo);
			DestroyMode(ma->L[ih]);
			free(ma->L[ih]);
		}

		DestroyModeArray(ma);	

	}
    DestroyGeometry(geo);
}
