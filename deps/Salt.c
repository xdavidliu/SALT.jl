#include "headers.h"



Mode *ReadModes(Geometry geo, char **namesin, char **namesout, int Nm){


	Mode *ms;

	int i;
	for(i=0; i<Nm; i++){

		double D;
		Mode m = ModeRead(namesin[i], geo, &D);

		if(i==0) geo->D = D;
		else if(D != geo->D)
			MyError("The input modes should all be at the same pump strength!");		

		sprintf(m->name, "%s", namesout[i]);


		addArrayMode(&ms, i, m);
	}
	return ms;

}


void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, double *eps, double *fprof, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
)
{
    Geometry geo;
    geo = CreateGeometry(N, M, h, Npml, Nc, LowerPML, eps, fprof, wa, y);    
    
    if(Dmax == 0.0){


		int i, added;
		Mode *ms = Passive(&added, BCPeriod, bl, k, wreal, wimag, modenorm, nev, geo);
		for(i = 0; i<added; i++){

			if(added == 1)
				sprintf(ms[i]->name, "%s", modeout);
			else 
				sprintf(ms[i]->name, "%s%i", modeout, i);


			Write(ms[i], geo);
			DestroyMode(ms[i]);
		}
		free(ms);

    }else{

		Mode *ms = ReadModes(geo, namesin, namesout, Nm);

		// hack, so that SALT.jl can access Creeper using just a pointer to modes.
        Creeper(dD, Dmax, thresholdw_tol, ftol, ms, printnewton, Nm, geo);    

		int ih;

		for(ih=0; ih<Nm; ih++){
			Write(ms[ih], geo);
			DestroyMode(ms[ih]);
			free(ms[ih]);
		}


	}
    DestroyGeometry(geo);
}
