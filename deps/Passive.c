#include "headers.h"

Mode CreateMode(Geometry geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]){
	Mode m = (Mode) malloc(sizeof(struct Mode_s) );
	CreateVec(2*Nxyzc(geo)+2, &m->vpsi);

	m->lasing = 0;
	m->ifix = ifix_;
	m->BCPeriod = BCPeriod_;
	int i, j;
	for(i=0; i<3; i++) for(j=0; j<2; j++) m->b[i][j] = b_[i][j];
	for(i=0; i<3; i++) m->k[i] = k_[i];

	return m;
} // make sure to change CopyMode too if this gets changed

void FillBop(Geometry geo, Mat Bop, dcomp w){
	VecCopy(geo->veps, geo->vscratch[0]);
	TimesI(geo, geo->vscratch[0], geo->vscratch[1]);

	Complexfun b;
	CreateComplexfun(&b, geo->vscratch[0], geo->vscratch[1]);
	int i; 
	for(i=b.ns; i<b.ne; i++)
		setc(&b, i, csqr(w)* valc(&b, i) );

	SetJacobian(geo, Bop, geo->vscratch[0], -2, 0, 0);
	SetJacobian(geo, Bop, geo->vscratch[1], -2, 1, 0);

	DestroyComplexfun(&b);

	AssembleMat(Bop);
}

// everything after modeout is directly from ReadGeometry
Mode *Passive(int *added, int BCPeriod, int *bl, double *k, double wreal, double wimag, double modenorm, int nev, Geometry geo){

    	tv t1, t2, t3;

	gettimeofday(&t1, NULL);

	int	i, b[3][2];

	for(i=0; i<3; i++){
		b[i][0] = bl[i];
		b[i][1] = 0;
	}

	Mat Mop;
	CreateSquareMatrix(2*Nxyzc(geo), 26, &Mop);

	MoperatorGeneralBlochFill(geo, Mop, b, BCPeriod, k, 0);
	AssembleMat(Mop);

	dcomp guess = -csqr(wreal + ComplexI * wimag);

	Mat Bop; CreateSquareMatrix(2*Nxyzc(geo), 2, &Bop);

	FillBop(geo, Bop, 1.0);

	EPS evps; 
        EPSCreate(PETSC_COMM_WORLD, &evps);
        EPSSetOperators(evps, Mop, Bop);
	EPSSetTarget(evps, creal(guess));
	EPSSetWhichEigenpairs(evps, EPS_TARGET_REAL);
	EPSSetDimensions(evps, nev,PETSC_DECIDE, PETSC_DECIDE);

        EPSSetFromOptions(evps);

	ST st; // doing this in case using Julia interface and cannot provide solver from command line arguments
	EPSGetST(evps, &st);
	STSetType(st, STSINVERT);

	KSP ksp;
	STGetKSP(st, &ksp);
	PC pc;
	KSPGetPC(ksp,&pc);
 	PCSetType(pc,PCLU);
  	PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);    

	gettimeofday(&t2, NULL);
        EPSSolve(evps);

        int nconv, j=0;
	EPSGetConverged(evps, &nconv);

	dcomp w; 

        Vec v, vi;
        MatGetVecs(Mop, &v, &vi);

	*added = 0;
	Mode *msp;

	if(nconv>0) for(i=0; i<nconv; i++){
		double lr, li;
		EPSGetEigenpair(evps, i, &lr, &li, v, vi);
		w = csqrt(-lr - ComplexI * li);

		//==============
		// tests if the eigenvector is conjugated. if unconjugated, vi should be
		// ( vI; -vR), i.e. the block version of -iv. Hence, if conjugated,
		// then || TimesI(v) + vi || > 1e-5 or something

		ScatterRange(v, geo->vscratch[0], 0, 0, xyzcrGrid(&geo->gN) );
		TimesI(geo, geo->vscratch[0], geo->vscratch[2]);

		VecSet(geo->vscratch[1], 0.0); // annoying last two elements
		ScatterRange(vi, geo->vscratch[1], 0, 0, xyzcrGrid(&geo->gN) );
		VecAXPY(geo->vscratch[2], -1.0, geo->vscratch[1]);

		double dvnorm, vnorm;
		VecNorm(geo->vscratch[2], NORM_2, &dvnorm);
		VecNorm(geo->vscratch[0], NORM_2, &vnorm);
		int fake = dvnorm / vnorm < 1e-5;
		//==============

		if( cimag(w) > 0.0 || fake ) continue; // could use ^ (XOR) here, but the positive frequency ones are exact duplicates.
		else j++;

		if(cimag(w) > 0.0) w = conj(w);

		Mode m = CreateMode(geo, 0, b, BCPeriod, k);
		ScatterRange(v, m->vpsi, 0, 0, xyzcrGrid(&geo->gN) );

		Fix(m, geo, modenorm);

		double psinorm, psifnorm;
		VecNorm(m->vpsi, NORM_2, &psinorm);
		VecPointwiseMult(geo->vscratch[0], m->vpsi, geo->vf);
		VecNorm(geo->vscratch[0], NORM_2, &psifnorm);

		PetscPrintf(PETSC_COMM_WORLD, "found mode #%i: w = %1.8g + (%1.8g) i; pumped fraction = %g\n", j-1, creal(w), cimag(w), psifnorm/psinorm );
		SetLast2(m->vpsi, creal(w), cimag(w) ); // make sure to get psinorm before doing this

		// TODO: free all Modes and Geometry since they are heap allocated when created

		addArrayMode(&msp, *added, m);
		(*added)++;

		if(nev == 1) break;
	}

	if(*added == 0) msp = NULL;

	EPSDestroy(&evps);
	MatDestroy(&Mop);
	MatDestroy(&Bop);
	VecDestroy(&v);	VecDestroy(&vi);

	gettimeofday(&t3, NULL);
	geo->D = 0.0; // to make sure pump strength for passive modes correct in Julia
	return msp;

//	PetscPrintf(PETSC_COMM_WORLD, "Formation and setup: %g s\nEPSSolve and Output: %g s\n", dt(t1, t2), dt(t2, t3) );
} 

