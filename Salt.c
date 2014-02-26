#include "headers.h"

#define MAXMODES 10
// TODO: allow more than 10 modes

PetscErrorCode ReadModes(ModeArray *ma, Geometry *geo, char namesin[MAXMODES][PETSC_MAX_PATH_LEN], char namesout[MAXMODES][PETSC_MAX_PATH_LEN], int Nm){


	int i;
	for(i=0; i<Nm; i++){

		double D;
		Mode* m = (Mode*) malloc(sizeof( Mode));
		ModeRead(m, namesin[i], geo, &D);

		if(i==0) geo->D = D;
		else if(D != geo->D)
			MyError("The input modes should all be at the same pump strength!");		

		sprintf(m->name, "%s", namesout[i]);

		if(i == 0) CreateModeArray(ma,m);
		else AddArrayMode(ma, m);

		Setup(m, geo);


	}
	return 0;

}




void FirstStep(ModeArray *mah, Mode *m, Geometry *geo, Vec vNh, Vec f, Vec dv, double c, double ftol, int printnewton){


	PetscPrintf(PETSC_COMM_WORLD, "Taking first step for mode \"%s\"...\n", m->name );

  int nh=0, ih;
  for(ih=0; ih<mah->size; ih++){ // find nh of m
		if( mah->L[ih] == m) break;
		else nh++;
  }


	if(vNh != m->vpsi){ // update vpsi's from v
		int ih =0;
		for(ih=0; ih<mah->size; ih++){
			ScatterRange((mah->L[ih])->vpsi, vNh, 0, ih*NJ(geo), NJ(geo) );
			
		}
	}

  while(1){

	if( LastProcess() ){ // try new c
		VecSetValue(vNh, offset(geo, nh)+Nxyzcr(geo)+1, c, INSERT_VALUES);	
		if( vNh != m->vpsi) VecSetValue(m->vpsi, Nxyzcr(geo)+1, c, INSERT_VALUES);
	}
	AssembleVec(vNh);
	AssembleVec(m->vpsi);


	double fnorm = FormJf(mah, geo, vNh, f, ftol, printnewton);


	if(  fnorm < ftol) break;
	
	KSPSolve( m->ksp, f, dv);
	if(printnewton) PetscPrintf(PETSC_COMM_WORLD, "\n");

	double dc = -GetValue(dv, offset(geo, nh)+Nxyzcr(geo)+1 );

	if( cabs(dc)/c < 0.5){
		VecAXPY(vNh, -1.0, dv);

		if(vNh != m->vpsi){ // update vpsi's from v
			int ih =0;
			for(ih=0; ih<mah->size; ih++){
				ScatterRange(vNh, (mah->L[ih])->vpsi, ih*NJ(geo), 0, NJ(geo) );
				
			}
		}
		// don't NewtonSolve here, that will be done immediately after
		break;
	}else if(c + dc < 0) c *= 0.5;
	else c = 0.5*(c + c+dc);

  }
  
  	PetscPrintf(PETSC_COMM_WORLD, "First step for mode \"%s\" complete!\n", m->name );  	


}




void Bundle(ModeArray *ma, Geometry *geo){




	int i, Nh = ma->size, Nj = 2*Nxyzc(geo)+2;
	if(Nh < 2) MyError("Bundle function is only for multimode!");
	
	Mat J; KSP ksp;
	CreateSquareMatrix( Nh*Nj, 0, &J);
	AllocateJacobian(J, geo);
	
	AddPlaceholders(J, geo);


	KSPCreate(PETSC_COMM_WORLD,&ksp);
	PC pc;
	KSPGetPC(ksp,&pc);
 	PCSetType(pc,PCLU);
  	PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
	// don't forget to change this in Setup too.

	KSPSetFromOptions(ksp);
	KSPSetOperators(ksp, J, J, SAME_PRECONDITIONER);
	// TODO: will probably want to merge all of this in with a generalized
	// multimode version of Mode::Setup
	

	for(i=0; i<SCRATCHNUM; i++){
		DestroyVec(&geo->vNhscratch[i]);
		MatGetVecs(J, &geo->vNhscratch[i], NULL);
	}


	int ih = 0;

	for(ih=0; ih<ma->size; ih++){
		Mode *m = ma->L[ih];
		
		DestroyMat( &m->J); // bundle shares J and v
		m->J = J;
		KSPDestroy(&m->ksp);
		m->ksp = ksp;
		
		MoperatorGeneralBlochFill(geo, J, m->b, m->BCPeriod, m->k, ih);
		AddRowDerivatives(J, geo, m->ifix, ih);
	}	
	
	AssembleMat(J);
	MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(J); 	
	
	
}



int FindModeAtThreshold(ModeArray *ma){

	int n = -1, ih;
	
	for(ih = 0; ih<ma->size; ih++){
		Mode *m = ma->L[ih];
		if( get_c(m) == 0.0 && m->lasing ){
			n = ih;
			break;
		}
	}
	return n;
	
}


// everything after Nm copied directly from ReadMode
void Salt(double dD, double Dmax, double thresholdw_tol, double ftol, char namesin[MAXMODES][PETSC_MAX_PATH_LEN], char namesout[MAXMODES][PETSC_MAX_PATH_LEN], int printnewton, int Nm, int N[3], int M[3], double h[3], int Npml[3], int Nc, int LowerPML, char *epsfile, char *fproffile, double wa, double y){

	Geometry Geo, *geo = &Geo;
	CreateGeometry(geo, N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);
	


	
  	ModeArray Ma, *ma = &Ma;




	// i is now number of modes read
	ReadModes(ma, geo, namesin, namesout, Nm);
	
	
    Vec f, dv;
    MatGetVecs( ma->L[0]->J, &dv, &f);




	int ih;
	for(; geo->D <= Dmax; geo->D = (geo->D+dD < Dmax? geo->D+dD: Dmax)){


	  	ModeArray Mah, *mah = &Mah;
		CreateFilter(ma, mah, 1); // lasing sub-array

	  
	  Vec vNh = ma->L[0]->vpsi, fNh = f, dvNh = dv;


	  if( mah->size > 0){ // lasing modes
	  
	
	  
	  	  int nt = FindModeAtThreshold(mah);
	  
		  if( nt != -1 && mah->size > 1){

		  	Bundle(mah, geo);
		  }

		  if(mah->size > 1){	 // these vectors will have been properly created in the last block
			vNh = geo->vNhscratch[2];
			fNh = geo->vNhscratch[3];
			dvNh = geo->vNhscratch[4];
		  }

		  if( nt != -1 ){
		  		geo->D += 0.5*dD;
				if(geo->D > Dmax) geo->D = Dmax;
				
		  		FirstStep(mah, mah->L[nt], geo, vNh, fNh, dvNh, 1.0, ftol, printnewton);
		  }
		  
		  
		  NewtonSolve(mah, geo,  vNh, fNh, dvNh, ftol, printnewton);  

		  
	  }

	  for(ih=0; ih<ma->size; ih++){ // now nonlasing modes
		Mode *m = ma->L[ih];
		if(m->lasing) continue;

		double wi_old = cimag(get_w(m));
		
		 ModeArray Ma_single, *ma_single = &Ma_single;
		 CreateModeArray(ma_single , m);
		
		NewtonSolve(ma_single , geo,  m->vpsi, f, dv, ftol, printnewton);
	  	DestroyModeArray(ma_single);
	  	
		double wi_new = cimag(get_w(m));

		if(wi_new > -thresholdw_tol && !m->lasing){
		
		
			ThresholdSearch(  wi_old, wi_new, geo->D-dD, geo->D, 
			mah, vNh, m, geo, f, dv, thresholdw_tol, ftol, printnewton); // todo: replace with vNh
			
		}
	  }

	  if(mah->size>0) DestroyModeArray(mah);	  
	  if(geo->D==Dmax) break;
	}

	for(ih=0; ih<ma->size; ih++){
		Write(ma->L[ih], geo);
		DestroyMode(ma->L[ih]);
		free(ma->L[ih]);
	}
	
  	DestroyModeArray(ma);	

	DestroyVec(&f);
	DestroyVec(&dv);
	DestroyGeometry(geo);


}


int main(int argc, char** argv){ 
	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); 


	double dD, Dmax, thresholdw_tol = OptionsDouble("-thresholdw_tol");	
	OptionsGetDouble("-dD", &dD);
	OptionsGetDouble("-Dmax", &Dmax);
	
	double ftol = OptionsDouble("-newtonf_tol");


	char optionin[PETSC_MAX_PATH_LEN] = "-in0",
	     optionout[PETSC_MAX_PATH_LEN] = "-out0",
		 namesin[MAXMODES][PETSC_MAX_PATH_LEN],
   		 namesout[MAXMODES][PETSC_MAX_PATH_LEN]; 
		

	int i=0;
	while(1){

		if( !OptionsGetString(optionin, namesin[i]) ) break;

		if( !OptionsGetString(optionout, namesout[i]) )
			MyError("number of -out less than number of -in!");
		
		i++;
		if(i > MAXMODES) MyError("exceeded mode limit!");
		sprintf(optionin, "-in%i", i);
		sprintf(optionout, "-out%i", i);
	}
	int Nm = i;

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


	int printnewton = OptionsInt("-printnewton");

	Salt(dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm,  N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);

	PetscPrintf(PETSC_COMM_WORLD, "\n");
	PetscPrintf(PETSC_COMM_WORLD, "TODO: a whole bunch of TODOs in Salt.c related to first step of multimode\n");	
	PetscPrintf(PETSC_COMM_WORLD, "future todo: add artificial crashes to enforce all the assumptions I'm making. For example, crash if any file read fails.\n");		
	PetscPrintf(PETSC_COMM_WORLD, "future todo: Make sure all MyError crashes crash all processes the way CHKERRQ does\n");		
	PetscPrintf(PETSC_COMM_WORLD, "TODO: compare multimode calculation speed to single mode with twice the pixels\n");		

	SlepcFinalize();	

}
