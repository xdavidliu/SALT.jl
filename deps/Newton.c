#include "salt.h"

// NOTE! as of 6/14, c has been redefined such that the hole burning term is 
// 1 / (1 + c^2 |psi|^2), i.e. all the rest of the garbage has been removed.
// Newton still works perfectly, because this amounts to a simple rescaling. For
// printing of data, we use the old convention. However, for outputting
// of c's to files, we're using the new convention. 

double EdgeIntensity(Mode m, Geometry geo){
// only works for sequential

	if( !m->lasing) return 0;

	TimesI(geo, m->vpsi, geo->vscratch[1]);
	Complexfun psi; Vecfun fprof;
	CreateComplexfun(&psi, m->vpsi, geo->vscratch[1]);
	CreateVecfun(&fprof, geo->vf);

	int i;
	double lastPsiSq = 0.0;
	for(i=0; i<Nxyzc(geo); i++){

		if( valr(&fprof, i) != 0.0 )
			lastPsiSq = sqr( cabs( valc(&psi, i) ) );
	} 

	DestroyComplexfun(&psi);
	DestroyVecfun(&fprof);

	return 2*sqr( get_c(m) ) * lastPsiSq;
}

void NewtonSolve(Mode *ms, Geometry geo, Vec v, Vec f, Vec dv, double ftol, int printnewton){
// f and dv are essentially scratch vectors.
// for L.size > 1, v is also essentially a scratch vector

	int ih =0, Nm;
	VecGetSize(v, &Nm); Nm /= NJ(geo);
	if( ms[0]->vpsi != v){  // update v from L's vpsi

		for(ih=0; ih<Nm; ih++){
			ScatterRange(ms[ih]->vpsi, v, 0, ih*NJ(geo), NJ(geo) );
		}
	}

	int its =0;
	tv t1, t2, t3;
	KSP ksp = ms[0]->ksp;

	VecSet(dv, 0.0);
	while(1){

		// removed stability hack for simplicity
		VecAXPY(v, -1.0, dv);

		if( ms[0]->vpsi != v && its !=0){ // update vpsi's from v
			for(ih=0; ih<Nm; ih++){
				ScatterRange(v, ms[ih]->vpsi, ih*NJ(geo), 0, NJ(geo) );
			}
		}

		gettimeofday(&t1, NULL);

		double fnorm = FormJf(ms, geo, v, f, ftol, printnewton);

		if(  fnorm < ftol)	break;

		gettimeofday(&t2, NULL);
		KSPSolve(ksp, f, dv);
		KSPSetOperators(ms[0]->ksp, ms[0]->J, ms[0]->J);

		gettimeofday(&t3, NULL);

		int printtime = 0; // disable printtime for now
		if(printtime)
		PetscPrintf(PETSC_COMM_WORLD, "formJ in %gs, solve in %gs\n", dt(t1, t2), dt(t2, t3) );
		else if(printnewton) PetscPrintf(PETSC_COMM_WORLD, "\n");

		its++;
		if(its > 6)
			KSPSetOperators(ms[0]->ksp, ms[0]->J, ms[0]->J);
		// refresh LU factorization. In next KSP solve, reset to SAME_PRECONDITIONER, do not need to worry about any other instances of KSPSolve
	}

	gettimeofday(&t2, NULL);

	if(printnewton){ 
		PetscPrintf(PETSC_COMM_WORLD, "\nconverged!\n" );
		for(ih=0; ih<Nm; ih++){
			dcomp w = get_w(ms[ih]);
			PetscPrintf(PETSC_COMM_WORLD, "%s at D = %g: w = %1.10g", ms[ih]->name,  geo->D, creal(w));
			if(! ms[ih]->lasing)
				PetscPrintf(PETSC_COMM_WORLD, " + i(%g)", cimag(w) );
			else
				PetscPrintf(PETSC_COMM_WORLD, ", c = %g", get_c(ms[ih]) / cabs(gamma_w(ms[ih], geo)) );
			// c in code redefined to absorb all the garbage multiplying the
			// |E|^2 in the hole burning. This outputted c uses the older
			// convention, on the other hand. -DL 6/14

			if( GetSize() == 1 && ms[ih]->lasing && geo->LowerPML==0 
			&& geo->gN.N[1] == 1 && geo->gN.N[2] == 1 && geo->Nc == 1)  
				PetscPrintf(PETSC_COMM_WORLD, ", |psi|^2_edge = %g", EdgeIntensity(ms[ih], geo) / sqr(cabs(gamma_w(ms[ih], geo))) );

			PetscPrintf(PETSC_COMM_WORLD, "\n");
		}
	}
}

void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, Mode *msh, Vec vNh, Vec fNh, Vec dvNh, Mode m, Geometry geo, Vec f, Vec dv, double ftol, int printnewton){

	dcomp mw = get_w(m);
	SetLast2(m->vpsi, creal(mw), 0.0);
	if(msh){ // lasing mode array passed in
		int N;
		VecGetSize(vNh, &N); 
	}

		
	// Threshold found if nonlasing residual with omega real is < ftol

	// 9/17/14, not sure why I had Nlasing > 1?vNh : m->vpsi for third argument here.
	// this is nonlasing mode, so don't care what Nlasing is.

	if( FormJf(&m, geo, m->vpsi, f, ftol, printnewton) < ftol ){
		PetscPrintf(PETSC_COMM_WORLD, "\nThreshold found for mode \"%s\" at D = %1.10g\n", m->name, geo->D);
		m->lasing = 1;
		return;
	}else{
		SetLast2(m->vpsi, creal(mw), cimag(mw));
	}


	geo->D = D_lo - (D_hi - D_lo)/(wimag_hi - wimag_lo) * wimag_lo;


	// set msh = NULL before calling ThresholdSearch if no lasing
	if( msh ) NewtonSolve(msh, geo, vNh, fNh, dvNh, ftol, printnewton);

	NewtonSolve(&m, geo, m->vpsi, f, dv, ftol, printnewton);
	mw = get_w(m);

	if( wimag_lo*wimag_hi > 0){ // both on same side of threshold
		if(wimag_lo > 0){ wimag_hi = wimag_lo; wimag_lo = cimag(mw); D_hi = D_lo; D_lo = geo->D;}
		else { wimag_lo = wimag_hi; wimag_hi = cimag(mw); D_lo = D_hi; D_hi = geo->D;}
	}else{// straddling threshold
		if(cimag(mw) > 0) {wimag_hi = cimag(mw); D_hi = geo->D;}
		else { wimag_lo = cimag(mw); D_lo = geo->D;}
	}

	if(printnewton)
	PetscPrintf(PETSC_COMM_WORLD, 
		"Searching... D=%g --> Im[w] = %g\n", geo->D, cimag(mw));
	ThresholdSearch(wimag_lo, wimag_hi, D_lo, D_hi, msh, vNh, fNh, dvNh, m, geo, f, dv, ftol, printnewton); 
}
