#include "salt.h"

int CreateFilter(Mode *ms, int size, int lasing, Mode **msp){
	int i, added =0;
	for(i=0; i<size; i++)
		if( ms[i]->lasing == lasing){
			addArrayMode(msp, added, ms[i]);
			added++;
		}

	if(!added) *msp = NULL; // so thresholdsearch receives correct argument
	return added;
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

PetscErrorCode Bundle(Mode *ms, int size, Geometry geo){
	int i, ih, Nh = size, Nj = 2*Nxyzc(geo)+2;
	if(Nh < 2) MyError("Bundle function is only for multimode!");

	Mat J; KSP ksp;
	CreateSquareMatrix( Nh*Nj, 0, &J);
	AllocateJacobian(J, geo);

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	PC pc;
	KSPGetPC(ksp,&pc);
 	PCSetType(pc,PCLU);
  	PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
	// don't forget to change this in Setup too.

	KSPSetFromOptions(ksp);
	KSPSetOperators(ksp, J, J); 
	// for petsc 3.4 and before, put SAME_PRECONDITIONER as fourth argument
	// for 3.5 and above, have only 3 arguments

	// TODO: will probably want to merge all of this in with a generalized
	// multimode version of Mode::Setup

	for(i=0; i<SCRATCHNUM; i++){
		if(geo->vNhscratch[i])
			VecDestroy(&geo->vNhscratch[i]);
		MatGetVecs(J, &geo->vNhscratch[i], NULL);
	}

	// when Bundle is called in the middle of Creeper, the modes
	// will have separate J's and ksps, so need to destroy each of them
	if(size == 2){ 
		for(ih=0; ih<size;ih++){
			Mode m = ms[ih];
			if(m->J){
				PetscErrorCode ierr = MatDestroy( &m->J); CHKERRQ(ierr);
			}
			if(m->ksp){
				KSPDestroy( &m->ksp);	
			}	
		}
	}else{ // size > 2; they will all have shared J and ksp
		   // except the one that just hit threshold
		   // but I'm destroying that one's J and ksp after
		   // thresholdSearch
		if( ms[0]->J ){
			PetscErrorCode ierr = MatDestroy( &ms[0]->J); CHKERRQ(ierr);
		}
		if( ms[0]->ksp) KSPDestroy( &ms[0]->ksp);
	}

	for(ih=0; ih<size; ih++){
		Mode m = ms[ih];
		m->J = J;// bundle shares J and v
		m->ksp = ksp;
		MoperatorGeneralBlochFill(geo, J, m->b, m->k, ih);
		AddRowDerivatives(J, geo, m->ifix, ih);
	}

	AssembleMat(J);
	MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(J); 
	return 0;
}

int FindModeAtThreshold(Mode *ms, int size){
	int n = -1, ih;

	for(ih = 0; ih<size; ih++){
		Mode m = ms[ih];
		if( get_c(m) == 0.0 && m->lasing ){
	
			if(n == -1)
				n = ih;
			else
				MyError("found two modes at threshold! FirstStep does not support multimode yet");
		}
	}
	return n;
}

void ComplexScale( Vec w, dcomp a, Vec scratch, Geometry geo){

	TimesI(geo, w, scratch);
	VecScale(w, creal(a));
	VecAXPY(w, cimag(a), scratch);

}

void ComplexPointwiseMult(Vec w, Vec u, Vec v, Vec scratch0, Vec scratch1, Geometry geo){
// w = u .* v
// i.e. wR = uR vR - uI vI
// wI = uR vI + uI vR

	int Nxyzc = xyzcGrid(&geo->gN);
	VecCopy(u, scratch0);
	ScatterRange(u, scratch0, 0, Nxyzc, Nxyzc );
	VecPointwiseMult(scratch0, scratch0, v);
	// scratch0 is now [uR vR; uR vI]
	
	VecCopy(u, scratch1);
	ScatterRange(u, scratch1, Nxyzc, 0, Nxyzc );
	// scratch1 is now [uI; uI]
	
	TimesI(geo, v, w);
	VecPointwiseMult(w, w, scratch1);
	// w is now [ -uI vI; uI vR]
	VecAXPY(w, 1.0, scratch0); 

}

// everything after Nm copied directly from ReadMode
int Creeper(double dD, double Dmax, double ftol, Mode *ms, int printnewton, int Nm, Geometry geo){
	
	double hugeval = 1.0e20;
	if(Dmax < 0.0) Dmax = hugeval; // hack, to instruct Creeper to stop upon threshold
	Mode *msh; // lasing mode subarray
	int ih, i, Nlasing, NlasingOld;

	Nlasing = CreateFilter(ms, Nm, 1, &msh);
	NlasingOld = Nlasing;
	if(Nlasing > 1) Bundle(msh, Nlasing, geo);

	for(ih=0; ih<Nm; ih++){
		if( !ms[ih]->lasing || Nlasing == 1)
			Setup( ms[ih], geo);
	}
    Vec f, dv;
//    MatGetVecs( ms[0]->J, &dv, &f); 
// 9/17/14: replaced MatGetVecs here with VecDuplicate. 2-mode lasing with 1-mode nonlasing no longer broken. Not sure why I didn't do this before.
	VecDuplicate(geo->vscratch[0], &dv);
	VecDuplicate(geo->vscratch[0], &f);

	for(; geo->D <= Dmax; geo->D = (geo->D+dD < Dmax? geo->D+dD: Dmax)){
	  	Nlasing = CreateFilter(ms, Nm, 1, &msh); // lasing sub-array
	  	Vec vNh = ms[0]->vpsi, fNh = f, dvNh = dv;
	 	if( Nlasing > 0){ // lasing modes  
	  	  int nt = FindModeAtThreshold(msh, Nlasing);
	  
		  // this is not called if we start from a threshold
		  if( nt != -1 && Nlasing > 1 && !msh[nt]->J ){
		  	Bundle(msh, Nlasing, geo);
		  }

		  if(Nlasing > 1){	 // these vectors will have been properly created in the last block
			vNh = geo->vNhscratch[2];
			fNh = geo->vNhscratch[3];
			dvNh = geo->vNhscratch[4];
		  }

		  if(Nlasing == 1)
		  	vNh = msh[0]->vpsi;

		  if( nt != -1 ){
		  		geo->D += 0.5*dD;
				if(geo->D > Dmax) geo->D = Dmax;
		  		FirstStep(msh, msh[nt], geo, vNh, fNh, dvNh, 1.0, ftol, printnewton);
		  }
		  NewtonSolve(msh, geo,  vNh, fNh, dvNh, ftol, printnewton);  
	  }

	
	  Mode mthreshold_nonlasing = 0;
	  for(ih=0; ih<Nm; ih++){ // now nonlasing modes

		Mode m = ms[ih];
		if(m->lasing || m == mthreshold_nonlasing) continue;
		double wi_old = cimag(get_w(m));

		NewtonSolve(&m, geo,  m->vpsi, f, dv, ftol, printnewton);

		double wi_new = cimag(get_w(m));

		if(wi_new > 0.0 && !m->lasing){

			ThresholdSearch(  wi_old, wi_new, geo->D-dD, geo->D, 
			msh, vNh, fNh, dvNh, m, geo, f, dv, ftol, printnewton);
			ih = -1; // reset to recalculate the rest of the lasing modes
			mthreshold_nonlasing = m;

			// now there will 2 or more lasing modes, so this threshold mode will join a multimode bundle
			if(Nlasing > 0){ 
				MatDestroy( &m->J);
				m->J = 0;
				KSPDestroy( &m->ksp);
				m->ksp = 0;
			}
		}
	  }

	  if(Nlasing>0) free(msh);	  
	  if(geo->D==Dmax || (Dmax == hugeval && mthreshold_nonlasing != 0) ) break;
	}

	// 070315: hack: output deps using perturbation theory and quadratic
	// programming method
	/*
	int output_deps = 0;
	PetscOptionsGetInt(PETSC_NULL,"-output_deps", &output_deps,NULL);
	if( output_deps == 1){
		PetscPrintf(PETSC_COMM_WORLD, "DEBUG: output_deps called!\n");

		ComplexPointwiseMult(geo->vscratch[0], ms[0]->vpsi, ms[1]->vpsi, geo->vscratch[1], geo->vscratch[2], geo);

		Output(geo->vscratch[0], "vpsi0psi1", "psi0psi1");

		VecCopy( ms[0]->vpsi, geo->vscratch[0]);
		ComplexScale( geo->vscratch[0], 2.0 - 3.0*ComplexI, geo->vscratch[1], geo);
		Output(geo->vscratch[0], "vpsi0a", "psi0a");

		Output(ms[0]->vpsi, "vpsi0", "psi0");
		Output(ms[1]->vpsi, "vpsi1", "psi1");

		// last two elements will be nonsensical, but no worries
	}
	*/
	//=======================

	PetscPrintf(PETSC_COMM_WORLD, "DEBUG: outputting H vector\n");
	Output(geo->vH, "Hvec", "H");

	PetscPrintf(PETSC_COMM_WORLD, "DEBUG: outputting eps vector\n");
	Output(geo->veps, "Epsvec", "Eps");	

	VecDestroy(&f);
	VecDestroy(&dv);

	Nlasing = CreateFilter(ms, Nm, 1, &msh);
	if(Nlasing > 0){
		MatDestroy(&msh[0]->J);
		KSPDestroy(&msh[0]->ksp);
		// these share the same J and ksp, so only one need to be destroyed
	}

	for(i=0; i<Nm; i++){
		if(!ms[i]->lasing || Nlasing == 1){
			MatDestroy(&ms[i]->J);
			KSPDestroy(&ms[i]->ksp);		
		}

		ms[i]->J = 0;
		ms[i]->ksp = 0;
	}

	for(i=0; i<SCRATCHNUM; i++){ // cleanup
		VecSet( geo->vH, 1.0);
		VecSet( geo->vscratch[i], 0.0);
		VecSet( geo->vMscratch[i], 0.0);
		if(geo->vNhscratch[i]){
			VecDestroy(&geo->vNhscratch[i]);
			geo->vNhscratch[i] = 0;
		}
	}
	return Nlasing - NlasingOld;
}
