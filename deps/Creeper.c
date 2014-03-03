#include "headers.h"

PetscErrorCode ReadModes(ModeArray ma, Geometry geo, char **namesin, char **namesout, int Nm){

	int i;
	for(i=0; i<Nm; i++){

		double D;
		Mode m = ModeRead(namesin[i], geo, &D);

		if(i==0) geo->D = D;
		else if(D != geo->D)
			MyError("The input modes should all be at the same pump strength!");		

		sprintf(m->name, "%s", namesout[i]);


		addArrayMode(&ma->L, ma->size, m);
		ma->size++;


	}
	return 0;

}

void FirstStep(Mode *ms, Mode m, Geometry geo, Vec vNh, Vec f, Vec dv, double c, double ftol, int printnewton){

	PetscPrintf(PETSC_COMM_WORLD, "Taking first step for mode \"%s\"...\n", m->name );

	int nh=0, ih, Nm;
	VecGetSize(vNh, &Nm); Nm /= NJ(geo);

	for(ih=0; ih<Nm; ih++){ // find nh of m
		if( ms[ih] == m) break;
		else nh++;
	}

	if(vNh != m->vpsi){ // update vpsi's from v
		int ih =0;
		for(ih=0; ih<Nm; ih++){
			ScatterRange((ms[ih])->vpsi, vNh, 0, ih*NJ(geo), NJ(geo) );
		
		}
	}

	while(1){

	if( LastProcess() ){ // try new c
		VecSetValue(vNh, offset(geo, nh)+Nxyzcr(geo)+1, c, INSERT_VALUES);	
		if( vNh != m->vpsi) VecSetValue(m->vpsi, Nxyzcr(geo)+1, c, INSERT_VALUES);
	}
	AssembleVec(vNh);
	AssembleVec(m->vpsi);

	double fnorm = FormJf(ms, geo, vNh, f, ftol, printnewton);

	if(  fnorm < ftol) break;

	KSPSolve( m->ksp, f, dv);
	if(printnewton) PetscPrintf(PETSC_COMM_WORLD, "\n");

	double dc = -GetValue(dv, offset(geo, nh)+Nxyzcr(geo)+1 );

	if( cabs(dc)/c < 0.5){
		VecAXPY(vNh, -1.0, dv);

		if(vNh != m->vpsi){ // update vpsi's from v
			int ih =0;
			for(ih=0; ih<Nm; ih++){
				ScatterRange(vNh, (ms[ih])->vpsi, ih*NJ(geo), 0, NJ(geo) );
			
			}
		}
		// don't NewtonSolve here, that will be done immediately after
		break;
	}else if(c + dc < 0) c *= 0.5;
	else c = 0.5*(c + c+dc);

	}

	PetscPrintf(PETSC_COMM_WORLD, "First step for mode \"%s\" complete!\n", m->name );  	

}

void Bundle(Mode *ms, int size, Geometry geo){

	int i, Nh = size, Nj = 2*Nxyzc(geo)+2;
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

	for(ih=0; ih<size; ih++){
		Mode m = ms[ih];
		
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

int FindModeAtThreshold(Mode *ms, int size){

	int n = -1, ih;

	for(ih = 0; ih<size; ih++){
		Mode m = ms[ih];
		if( get_c(m) == 0.0 && m->lasing ){
			n = ih;
			break;
		}
	}
	return n;

}

// everything after Nm copied directly from ReadMode
void Creeper(double dD, double Dmax, double thresholdw_tol, double ftol, Mode *ms, int printnewton, int Nm, Geometry geo){


	ModeArray ma = CreateModeArray();
	int ih;
	for(ih=0; ih<Nm; ih++){

		addArrayMode(&ma->L, ma->size, ms[ih]);
		ma->size++;

		Setup( ms[ih], geo); // TODO: bundle if multiple lasing modes
	}


    Vec f, dv;
    MatGetVecs( ma->L[0]->J, &dv, &f);


	for(; geo->D <= Dmax; geo->D = (geo->D+dD < Dmax? geo->D+dD: Dmax)){

	  	ModeArray mah =	CreateFilter(ma, 1); // lasing sub-array

	  
	  Vec vNh = ma->L[0]->vpsi, fNh = f, dvNh = dv;

	  if( mah->size > 0){ // lasing modes
	  

	  
	  	  int nt = FindModeAtThreshold(mah->L, mah->size);
	  
		  if( nt != -1 && mah->size > 1){

		  	Bundle(mah->L, mah->size, geo);
		  }

		  if(mah->size > 1){	 // these vectors will have been properly created in the last block
			vNh = geo->vNhscratch[2];
			fNh = geo->vNhscratch[3];
			dvNh = geo->vNhscratch[4];
		  }

		  if(mah->size == 1)
		  	vNh = mah->L[0]->vpsi;

		  if( nt != -1 ){
		  		geo->D += 0.5*dD;
				if(geo->D > Dmax) geo->D = Dmax;
				
		  		FirstStep(mah->L, mah->L[nt], geo, vNh, fNh, dvNh, 1.0, ftol, printnewton);
		  }
		  
		  
		  NewtonSolve(mah->L, geo,  vNh, fNh, dvNh, ftol, printnewton);  

		  
	  }

	  for(ih=0; ih<ma->size; ih++){ // now nonlasing modes
		Mode m = ma->L[ih];
		if(m->lasing) continue;

		double wi_old = cimag(get_w(m));
		
		 ModeArray ma_single = CreateModeArray();
		addArrayMode(&ma_single->L, 0, m);
		ma_single->size = 1;

	

		
		NewtonSolve(ma_single->L, geo,  m->vpsi, f, dv, ftol, printnewton);
	  	DestroyModeArray(ma_single);
	  	
		double wi_new = cimag(get_w(m));

		if(wi_new > -thresholdw_tol && !m->lasing){
		
			// set msh to NULL if no lasing modes
			Mode *msh = ( mah->size > 0? mah->L : NULL);

			ThresholdSearch(  wi_old, wi_new, geo->D-dD, geo->D, 
			msh, vNh, m, geo, f, dv, thresholdw_tol, ftol, printnewton);
			
		}
	  }

	  if(mah->size>0) DestroyModeArray(mah);	  
	  if(geo->D==Dmax) break;
	}



	DestroyVec(&f);
	DestroyVec(&dv);

}