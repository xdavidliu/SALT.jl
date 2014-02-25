#include "headers.h"


double EdgeIntensity(Mode &m, Geometry *geo){

	if(geo->gN.N[1] != 1 || geo->gN.N[2] != 1 || geo->LowerPML != 0 || geo->Nc != 1)
		MyError("EdgeIntensity is only for 1d symmetric TM fields!");

	if( !m.lasing) return 0;
	
	double psiR = GetValue(m.vpsi, Mxyz(geo)-1 ),
			psiI = GetValue(m.vpsi, Mxyz(geo)-1+Nxyzc(geo) );

	return 2*sqr( getc(&m) ) * (  sqr(psiR ) + sqr(psiI) );

}



void NewtonSolve(modelist &L, Geometry *geo, Vec v, Vec f, Vec dv){
// f and dv are essentially scratch vectors.
// for L.size > 1, v is also essentially a scratch vector

	if( (*L.begin())->vpsi != v){  // update v from L's vpsi
		int ih =0;
		FORMODES(L, it){
			ScatterRange((*it)->vpsi, v, 0, ih*NJ(geo), NJ(geo) );
			ih++;
		}
	}


	int its =0;
	tv t1, t2, t3;
	KSP ksp = (*L.begin() )->ksp;

	VecSet(dv, 0.0);
	while(true){

		// removed stability hack for simplicity	
		VecAXPY(v, -1.0, dv);	

		if( (*L.begin())->vpsi != v && its !=0){ // update vpsi's from v
			int ih =0;
			FORMODES(L, it){
				ScatterRange(v, (*it)->vpsi, ih*NJ(geo), 0, NJ(geo) );
				ih++;
			}
		}
		

		gettimeofday(&t1, NULL);
		if( FormJf(L, geo, v, f) < OptionsDouble("-newtonf_tol"))	break;
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
		FORMODES(L, it){
			dcomp w = getw(*it);
			PetscPrintf(PETSC_COMM_WORLD, "%s at D = %g: w = %g + i(%g)", 
				(*it)->name,  geo->D, w.real(), w.imag());
				
			if( (*it)->lasing )  PetscPrintf(PETSC_COMM_WORLD, ", |psi|^2_edge = %g", EdgeIntensity( **it, geo));
				
			PetscPrintf(PETSC_COMM_WORLD, "\n");	
		}
	}

}




void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, modelist &Lh, Vec vNh, Mode& m, Geometry *geo, Vec f, Vec dv){

		
	modelist L;
	L.push_back(&m);
	
	dcomp mw = getw(&m);
	if( std::abs(mw.imag()) < OptionsDouble("-thresholdw_tol") ){
		SetLast2(m.vpsi, mw.real(), 0.0);
		PetscPrintf(PETSC_COMM_WORLD, "Threshold found for mode \"%s\" at D = %1.10g\n", m.name, geo->D);
		m.lasing = 1;
		return;
	}

	geo->D = D_lo - (D_hi - D_lo)/(wimag_hi - wimag_lo) * wimag_lo;
	if(Lh.size() > 0) NewtonSolve(Lh, geo, vNh, f, dv);

		// if searching a single mode with no lasing, pass empty list
	NewtonSolve(L, geo, m.vpsi, f, dv);
	mw = getw(&m);
	
	if( wimag_lo*wimag_hi > 0){ // both on same side of threshold
		if(wimag_lo > 0){ wimag_hi = wimag_lo; wimag_lo = mw.imag(); D_hi = D_lo; D_lo = geo->D;}
		else { wimag_lo = wimag_hi; wimag_hi = mw.imag(); D_lo = D_hi; D_hi = geo->D;}
	}else{// straddling threshold
		if(mw.imag() > 0) {wimag_hi = mw.imag(); D_hi = geo->D;}
		else { wimag_lo = mw.imag(); D_lo = geo->D;}
	}


	if(OptionsInt("-printnewton"))
	PetscPrintf(PETSC_COMM_WORLD, 
		"Searching... D=%g --> Im[w] = %g\n", geo->D, mw.imag());
	ThresholdSearch(wimag_lo, wimag_hi, D_lo, D_hi, Lh, vNh, m, geo, f, dv); 	
}
