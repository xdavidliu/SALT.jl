#include "headers.h"




int main(int argc, char** argv){ 
	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); 


	double Dmax;
	OptionsGetDouble("-Dmax", &Dmax);
	


	// ======== copied directly from ReadGeometry ======== //

	int N[3], M[3], Npml[3], Nc, LowerPML;
	double h[3];

	OptionsXYZInt("-N", N);
	OptionsXYZInt("-M", M);

	OptionsXYZInt("-Npml", Npml);
	OptionsXYZDouble("-h", h);

	OptionsGetInt("-Nc", &Nc);
	OptionsGetInt("-LowerPML", &LowerPML);


	char epsfile[PETSC_MAX_PATH_LEN], fproffile[PETSC_MAX_PATH_LEN];

	OptionsGetString("-epsfile", epsfile);
	OptionsGetString("-fproffile", fproffile);

	double wa, y;
	OptionsGetDouble("-wa", &wa);
	OptionsGetDouble("-gamma", &y);

	// ======== copied directly from ReadGeometry ======== //

	Geometry Geo, *geo = &Geo;
	CreateGeometry(geo, N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);


	if(Dmax == 0.0){ // Passive
		double wguess_real, wguess_imag, modenorm = OptionsDouble("-norm");
		OptionsGetDouble("-wreal", &wguess_real);
		OptionsGetDouble("-wimag", &wguess_imag);
	
		int bl[3], BCPeriod, nev = OptionsInt("-nev");
		OptionsGetInt("-BCPeriod", &BCPeriod);
		OptionsXYZInt("-b", bl);
	
		double k[3] = {0, 0, 0};
		OptionsXYZDouble("-k", k);
	
		char s[PETSC_MAX_PATH_LEN];
		OptionsGetString("-passiveout", s);
		Passive(BCPeriod, bl, k, wguess_real, wguess_imag, modenorm, nev, s, geo);

	}else{ // Creeper


		double dD, thresholdw_tol = OptionsDouble("-thresholdw_tol");	
		OptionsGetDouble("-dD", &dD);
	
		double ftol = OptionsDouble("-newtonf_tol");


		char optionin[PETSC_MAX_PATH_LEN] = "-in0",
			 optionout[PETSC_MAX_PATH_LEN] = "-out0",
			 name[PETSC_MAX_PATH_LEN], name2[PETSC_MAX_PATH_LEN],
			 **namesin,
			 **namesout;
		

		int i=0;
		while(1){

			if( !OptionsGetString(optionin, name) ) break;
			if( !OptionsGetString(optionout, name2) )
				MyError("number of -out less than number of -in!");
		

if(i==0){
	namesin = (char **) malloc( 1*sizeof(char*) );
	namesout = (char **) malloc( 1*sizeof(char*) );
}else{
	namesin = (char **) realloc( namesin, (i+1)*sizeof(char*) );
	namesout = (char **) realloc( namesout, (i+1)*sizeof(char*) );
}

	namesin[i] = (char *) malloc( PETSC_MAX_PATH_LEN * sizeof(char) );
	namesout[i] = (char *) malloc( PETSC_MAX_PATH_LEN * sizeof(char) );

sprintf(namesin[i], "%s", name);
sprintf(namesout[i], "%s", name2);

	i++;

			sprintf(optionin, "-in%i", i);
			sprintf(optionout, "-out%i", i);
		}
		int Nm = i;

		int printnewton = OptionsInt("-printnewton");

		Creeper(dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm, geo);

		for(i=0; i<Nm; i++){
			free( namesin[i]);
			free( namesout[i]);
		}
		free(namesin);
		free(namesout);

	}


	DestroyGeometry(geo);

	PetscPrintf(PETSC_COMM_WORLD, "\n");
	PetscPrintf(PETSC_COMM_WORLD, "TODO: a whole bunch of TODOs in Salt.c related to first step of multimode\n");	
	PetscPrintf(PETSC_COMM_WORLD, "future todo: add artificial crashes to enforce all the assumptions I'm making. For example, crash if any file read fails.\n");		
	PetscPrintf(PETSC_COMM_WORLD, "future todo: Make sure all MyError crashes crash all processes the way CHKERRQ does\n");		
	PetscPrintf(PETSC_COMM_WORLD, "TODO: compare multimode calculation speed to single mode with twice the pixels\n");		

	SlepcFinalize();	

}
