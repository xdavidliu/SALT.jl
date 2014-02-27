#include <slepc.h>


int OptionsGetString(const char* c, char* a){ 
	PetscBool flg;
	PetscOptionsGetString(PETSC_NULL,c, a, PETSC_MAX_PATH_LEN, &flg); 
	return flg;
}

int OptionsGetInt(const char* c, int* a){
	PetscBool flg;
	PetscOptionsGetInt(PETSC_NULL,c,a,&flg);
	return flg;
}


void OptionsXYZInt(const char* prefix, int* a){
	int i;
	char option[PETSC_MAX_PATH_LEN];
	const char x[3] = {'x', 'y', 'z'};
	for(i=0; i<3; i++){
		sprintf(option, "%s%c", prefix, x[i]);
		OptionsGetInt(option, &a[i]);
	}

}

int OptionsGetDouble(const char* c, double* a){
	PetscBool flg;
	PetscOptionsGetReal(PETSC_NULL,c,a,&flg); 
	return flg;
}


void OptionsXYZDouble(const char* prefix, double* a){
	int i;
	char option[PETSC_MAX_PATH_LEN];
	const char x[3] = {'x', 'y', 'z'};
	for(i=0; i<3; i++){
		sprintf(option, "%s%c", prefix, x[i]);
		OptionsGetDouble(option, &a[i]);
	}

}

PetscErrorCode MyError(const char* message);
void CreateVec(int N, Vec *x);
void ReadVectorC(FILE *fp, int N, Vec v);

void Salt(int *N, int *M, double *h, int *Npml, int Nc, int LowerPML, double *eps, double *fprof, double wa, double y,  // <-- Geometry parameters
int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, char *modeout,  // <--- Passive parameters
double dD, double Dmax, double thresholdw_tol, double ftol, char **namesin, char **namesout, int printnewton, int Nm // <--- Creeper parameters
);


int main(int argc, char** argv){ 
	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); 


	double Dmax;
	OptionsGetDouble("-Dmax", &Dmax);
	


	// ======== copied directly from ReadGeometry ======== //

	int N[3], M[3], Npml[3], Nc, LowerPML, i;
	double h[3];

	OptionsXYZInt("-N", N);
	OptionsXYZInt("-M", M);

	OptionsXYZInt("-Npml", Npml);
	OptionsXYZDouble("-h", h);

	OptionsGetInt("-Nc", &Nc);
	OptionsGetInt("-LowerPML", &LowerPML);




	double wa, y;
	OptionsGetDouble("-wa", &wa);
	OptionsGetDouble("-gamma", &y);

	// ======== copied directly from ReadGeometry ======== //



	double wreal = 0., wimag=0., modenorm=1., 
		k[3] = {0}, dD = 0.0, thresholdw_tol=0., ftol = 0.0;
	int bl[3] = {0}, BCPeriod=0, nev=0, Nm = 0, printnewton = 0;
	char modeout[PETSC_MAX_PATH_LEN] = "", **namesin = NULL, **namesout = NULL;

	if(Dmax == 0.0){ // Passive
		OptionsGetDouble("-norm", &modenorm);
		OptionsGetDouble("-wreal", &wreal);
		OptionsGetDouble("-wimag", &wimag);
	

		OptionsGetInt("-nev", &nev);
		OptionsGetInt("-BCPeriod", &BCPeriod);
		OptionsXYZInt("-b", bl);
		OptionsXYZDouble("-k", k);
	
		OptionsGetString("-passiveout", modeout);

	}else{ // Creeper


		OptionsGetDouble("-thresholdw_tol", &thresholdw_tol);	
		OptionsGetDouble("-dD", &dD);
	
		OptionsGetDouble("-newtonf_tol", &ftol);


		char optionin[PETSC_MAX_PATH_LEN] = "-in0",
			 optionout[PETSC_MAX_PATH_LEN] = "-out0",
			 name[PETSC_MAX_PATH_LEN], name2[PETSC_MAX_PATH_LEN];
		

		i=0;
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
		Nm = i;

		OptionsGetInt("-printnewton", &printnewton);
	}

// ======== read eps and fprof from file ======= //


	char epsfile[PETSC_MAX_PATH_LEN], fproffile[PETSC_MAX_PATH_LEN];

	OptionsGetString("-epsfile", epsfile);
	OptionsGetString("-fproffile", fproffile);

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


Salt(N, M, h, Npml, Nc, LowerPML, eps, fprof, wa, y,
BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout,
dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm);

	VecRestoreArray(veps, &eps);
	VecRestoreArray(vfprof, &fprof);

	VecDestroy(&veps);
	VecDestroy(&vfprof);


	for(i=0; i<Nm; i++){
		free( namesin[i]);
		free( namesout[i]);
	}
	free(namesin);
	free(namesout);


/*
	PetscPrintf(PETSC_COMM_WORLD, "\n");
	PetscPrintf(PETSC_COMM_WORLD, "TODO: a whole bunch of TODOs in Salt.c related to first step of multimode\n");	
	PetscPrintf(PETSC_COMM_WORLD, "future todo: add artificial crashes to enforce all the assumptions I'm making. For example, crash if any file read fails.\n");		
	PetscPrintf(PETSC_COMM_WORLD, "future todo: Make sure all MyError crashes crash all processes the way CHKERRQ does\n");		
	PetscPrintf(PETSC_COMM_WORLD, "TODO: compare multimode calculation speed to single mode with twice the pixels\n");		
*/

	SlepcFinalize();	

}