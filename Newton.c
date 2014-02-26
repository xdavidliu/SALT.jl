#include "headers.h"


double EdgeIntensity(Mode *m, Geometry *geo){

	if(geo->gN.N[1] != 1 || geo->gN.N[2] != 1 || geo->LowerPML != 0 || geo->Nc != 1)
		MyError("EdgeIntensity is only for 1d symmetric TM fields!");

	if( !m->lasing) return 0;
	
	double psiR = GetValue(m->vpsi, Mxyz(geo)-1 ),
			psiI = GetValue(m->vpsi, Mxyz(geo)-1+Nxyzc(geo) );

	return 2*sqr( get_c(m) ) * (  sqr(psiR ) + sqr(psiI) );

}



void NewtonSolve(ModeArray *ma, Geometry *geo, Vec v, Vec f, Vec dv){
// f and dv are essentially scratch vectors.
// for L.size > 1, v is also essentially a scratch vector

	int ih =0;
	if( ma->L[0]->vpsi != v){  // update v from L's vpsi

		
		for(ih=0; ih<ma->size; ih++){
			ScatterRange(ma->L[ih]->vpsi, v, 0, ih*NJ(geo), NJ(geo) );
		}
	}


	int its =0;
	tv t1, t2, t3;
	KSP ksp = ma->L[0]->ksp;

	VecSet(dv, 0.0);
	while(1){

		// removed stability hack for simplicity	
		VecAXPY(v, -1.0, dv);	

		if( ma->L[0]->vpsi != v && its !=0){ // update vpsi's from v
			for(ih=0; ih<ma->size; ih++){
				ScatterRange(v, ma->L[ih]->vpsi, ih*NJ(geo), 0, NJ(geo) );
			}
		}
		

		gettimeofday(&t1, NULL);
		


		double fnorm = FormJf(ma, geo, v, f);
		
		if(  fnorm < OptionsDouble("-newtonf_tol"))	break;
		gettimeofday(&t2, NULL);
		KSPSolve(ksp, f, dv);
		gettimeofday(&t3, NULL);


		if(OptionsInt("-printtime"))
		PetscPrintf(PETSC_COMM_WORLD, "formJ in %gs, solve in %gs\n", dt(t1, t2), dt(t2, t3) );
		else PetscPrintf(PETSC_COMM_WORLD, "\n");

		its++;
		if(its > 8) MyError("Something's wrong, Newton shouldn't take this long to converge. Try a different starting point.\n");
	}

	gettimeofday(&t2, NULL);


	
	// removed print statement; redo these for multimode.
	if(OptionsInt("-printnewton")){ 
		PetscPrintf(PETSC_COMM_WORLD, "\nconverged!\n" );
		for(ih=0; ih<ma->size; ih++){
			dcomp w = get_w(ma->L[ih]);
			PetscPrintf(PETSC_COMM_WORLD, "%s at D = %g: w = %g + i(%g)", 
				ma->L[ih]->name,  geo->D, creal(w), cimag(w));
				
			if( ma->L[ih]->lasing && geo->LowerPML==0 )  PetscPrintf(PETSC_COMM_WORLD, ", |psi|^2_edge = %g", EdgeIntensity(ma->L[ih], geo));
				
			PetscPrintf(PETSC_COMM_WORLD, "\n");	
		}
	}

}




void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, ModeArray *mah, Vec vNh, Mode *m, Geometry *geo, Vec f, Vec dv){

	
	dcomp mw = get_w(m);
	if( cabs(cimag(mw)) < OptionsDouble("-thresholdw_tol") ){
		SetLast2(m->vpsi, creal(mw), 0.0);
		PetscPrintf(PETSC_COMM_WORLD, "Threshold found for mode \"%s\" at D = %1.10g\n", m->name, geo->D);
		m->lasing = 1;
		return;
	}



	geo->D = D_lo - (D_hi - D_lo)/(wimag_hi - wimag_lo) * wimag_lo;
	if(mah->size>0) NewtonSolve(mah, geo, vNh, f, dv);

	ModeArray Ma, *ma = &Ma;	
	CreateModeArray(ma, m);
	
		// if searching a single mode with no lasing, pass empty list
	NewtonSolve(ma, geo, m->vpsi, f, dv);
	DestroyModeArray(ma);

	mw = get_w(m);
	
	if( wimag_lo*wimag_hi > 0){ // both on same side of threshold
		if(wimag_lo > 0){ wimag_hi = wimag_lo; wimag_lo = cimag(mw); D_hi = D_lo; D_lo = geo->D;}
		else { wimag_lo = wimag_hi; wimag_hi = cimag(mw); D_lo = D_hi; D_hi = geo->D;}
	}else{// straddling threshold
		if(cimag(mw) > 0) {wimag_hi = cimag(mw); D_hi = geo->D;}
		else { wimag_lo = cimag(mw); D_lo = geo->D;}
	}


	if(OptionsInt("-printnewton"))
	PetscPrintf(PETSC_COMM_WORLD, 
		"Searching... D=%g --> Im[w] = %g\n", geo->D, cimag(mw));
	ThresholdSearch(wimag_lo, wimag_hi, D_lo, D_hi, mah, vNh, m, geo, f, dv); 	
}
