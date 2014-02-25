#include "headers.h"



void FillBop(Geometry& geo, Mat Bop, dcomp w){


	VecCopy(geo.veps, geo.vscratch[0]);
	geo.TimesI(geo.vscratch[0], geo.vscratch[1]);
	
	Complexfun b;
	CreateComplexfun(&b, geo.vscratch[0], geo.vscratch[1]);

	for(int i=b.ns; i<b.ne; i++)
		setc(&b, i, sqr(w)* valc(&b, i) );

	geo.SetJacobian(Bop, geo.vscratch[0], -2, 0, 0);
	geo.SetJacobian(Bop, geo.vscratch[1], -2, 1, 0);	
	
	DestroyComplexfun(&b);
	
	AssembleMat(Bop);

}


int main(int argc, char** argv){ SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); {



    	tv t1, t2, t3;
	Geometry geo;
	gettimeofday(&t1, NULL);

	int	b[3][2], BCPeriod, bl[3];
	OptionsGetInt("-BCPeriod", &BCPeriod);
	OptionsXYZInt("-b", bl);

	for(int i=0; i<3; i++){ b[i][0]=bl[i]; b[i][1] = 0;}

        Mat Mop;
	CreateSquareMatrix(2*Nxyzc(&geo), 26, &Mop);




	double k[3] = {0, 0, 0};
	OptionsXYZDouble("-k", k);




	MoperatorGeneralBlochFill(&geo, Mop, b, BCPeriod, k);
	AssembleMat(Mop);









	double wguess_real, wguess_imag;
	OptionsGetDouble("-wreal", &wguess_real);
	OptionsGetDouble("-wimag", &wguess_imag);
	dcomp guess = -sqr(wguess_real + ComplexI * wguess_imag);


	Mat Bop; CreateSquareMatrix(2*Nxyzc(&geo), 2, &Bop);



	FillBop(geo, Bop, 1.0);







	EPS evps; 
        EPSCreate(PETSC_COMM_WORLD, &evps);
        EPSSetOperators(evps, Mop, Bop);
	EPSSetTarget(evps, guess.real());
	EPSSetWhichEigenpairs(evps, EPS_TARGET_REAL);


	

        EPSSetFromOptions(evps);
    


	gettimeofday(&t2, NULL);
        EPSSolve(evps);


	


        int nconv, j=0;
	EPSGetConverged(evps, &nconv);
	
	dcomp w; 

        Vec v, vi;
        MatGetVecs(Mop, &v, &vi);

	if(nconv>0) for(int i=0; i<nconv; i++){
		double lr, li;
		EPSGetEigenpair(evps, i, &lr, &li, v, vi);
		w = std::sqrt(-lr - ComplexI * li);
		
		//==============
		// tests if the eigenvector is conjugated. if unconjugated, vi should be
		// ( vI; -vR), i.e. the block version of -iv. Hence, if conjugated,
		// then || TimesI(v) + vi || > 1e-5 or something
		

		ScatterRange(v, geo.vscratch[0], 0, 0, xyzcrGrid(&geo.gN) );
		geo.TimesI(geo.vscratch[0], geo.vscratch[2]);
		
		VecSet(geo.vscratch[1], 0.0); // annoying last two elements
		ScatterRange(vi, geo.vscratch[1], 0, 0, xyzcrGrid(&geo.gN) );
		VecAXPY(geo.vscratch[2], -1.0, geo.vscratch[1]);
		
		double dvnorm, vnorm;
		VecNorm(geo.vscratch[2], NORM_2, &dvnorm);
		VecNorm(geo.vscratch[0], NORM_2, &vnorm);
		int fake = dvnorm / vnorm < 1e-5;
		//==============

		
		if( w.imag() > 0.0 || fake ) continue; // could use ^ (XOR) here, but the positive frequency ones are exact duplicates.
		else j++;


		if(w.imag() > 0.0) w = std::conj(w);


		Mode m(geo, 0, b, BCPeriod, k);
		ScatterRange(v, m.vpsi, 0, 0, xyzcrGrid(&geo.gN) );


		Fix(&m, geo);
		
		double psinorm, psifnorm;
		VecNorm(m.vpsi, NORM_2, &psinorm);
		VecPointwiseMult(geo.vscratch[0], m.vpsi, geo.vf);
		VecNorm(geo.vscratch[0], NORM_2, &psifnorm);

		PetscPrintf(PETSC_COMM_WORLD, "found mode #%i: w = %1.8g + (%1.8g) i; pumped fraction = %g\n", j, w.real(), w.imag(), psifnorm/psinorm );
		SetLast2(m.vpsi, w.real(), w.imag() ); // make sure to get psinorm before doing this


		char s[PETSC_MAX_PATH_LEN];
		OptionsGetString("-passiveout", s);
		sprintf(m.name, "%s%i", s, j);
		m.Write(geo);
	}
	EPSDestroy(&evps);
	DestroyMat(&Mop);	DestroyMat(&Bop);
	DestroyVec(&v);	DestroyVec(&vi);


	gettimeofday(&t3, NULL);	
	
	PetscPrintf(PETSC_COMM_WORLD, "Formation and setup: %g s\nEPSSolve and Output: %g s\n", 
		dt(t1, t2), dt(t2, t3) );

		
} SlepcFinalize();	}
