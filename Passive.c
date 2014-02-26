#include "headers.h"



void FillBop(Geometry *geo, Mat Bop, dcomp w){


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
void Passive(int BCPeriod, int bl[3], double k[3], double wreal, double wimag, double modenorm, int nev, char *modeout, int N[3], int M[3], double h[3], int Npml[3], int Nc, int LowerPML, char *epsfile, char *fproffile, double wa, double y){

    	tv t1, t2, t3;
	Geometry Geo;
	Geometry *geo = &Geo; // to make consistent with other functions
	CreateGeometry(geo, N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);	

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


		Mode M, *m = &M;
		CreateMode(m, geo, 0, b, BCPeriod, k);
		ScatterRange(v, m->vpsi, 0, 0, xyzcrGrid(&geo->gN) );


		Fix(m, geo, modenorm);
		
		double psinorm, psifnorm;
		VecNorm(m->vpsi, NORM_2, &psinorm);
		VecPointwiseMult(geo->vscratch[0], m->vpsi, geo->vf);
		VecNorm(geo->vscratch[0], NORM_2, &psifnorm);

		PetscPrintf(PETSC_COMM_WORLD, "found mode #%i: w = %1.8g + (%1.8g) i; pumped fraction = %g\n", j, creal(w), cimag(w), psifnorm/psinorm );
		SetLast2(m->vpsi, creal(w), cimag(w) ); // make sure to get psinorm before doing this



		sprintf(m->name, "%s%i", modeout, j);
		Write(m, geo);
		DestroyMode(m);
	}
	EPSDestroy(&evps);
	DestroyMat(&Mop);	DestroyMat(&Bop);
	DestroyVec(&v);	DestroyVec(&vi);
	DestroyGeometry(geo);

	gettimeofday(&t3, NULL);	
	
	PetscPrintf(PETSC_COMM_WORLD, "Formation and setup: %g s\nEPSSolve and Output: %g s\n", 
		dt(t1, t2), dt(t2, t3) );

} 


int main(int argc, char** argv){ 

	SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); 


	double wguess_real, wguess_imag;
	OptionsGetDouble("-wreal", &wguess_real);
	OptionsGetDouble("-wimag", &wguess_imag);

	int bl[3], BCPeriod;
	OptionsGetInt("-BCPeriod", &BCPeriod);
	OptionsXYZInt("-b", bl);

	double k[3] = {0, 0, 0};
	OptionsXYZDouble("-k", k);

	char s[PETSC_MAX_PATH_LEN];
	OptionsGetString("-passiveout", s);

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

	double wa, y, modenorm = OptionsDouble("-norm");
	OptionsGetDouble("-wa", &wa);
	OptionsGetDouble("-gamma", &y);
	
	int nev = OptionsInt("-nev");

	// ======== copied directly from ReadGeometry ======== //

	Passive(BCPeriod, bl, k, wguess_real, wguess_imag, modenorm, nev, s,N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y);

		
SlepcFinalize();	
}
