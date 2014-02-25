#include "headers.h"




Mode::Mode(Geometry& geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]){

	CreateVec(2*Nxyzc(&geo)+2, &vpsi);
	CreateSquareMatrix(2*Nxyzc(&geo)+2, 0, &J);
	KSPCreate(PETSC_COMM_WORLD,&ksp);

	lasing = 0;
	ifix = ifix_;
	BCPeriod = BCPeriod_;
	for(int i=0; i<3; i++) for(int j=0; j<2; j++) b[i][j] = b_[i][j];
	for(int i=0; i<3; i++) k[i] = k_[i];
}



void ModeRead(Mode *m, char *Name, Geometry& geo, double *Dout){

	sprintf(m->name, "%s", Name );
	CreateVec(2*Nxyzc(&geo)+2, &m->vpsi);
	CreateSquareMatrix(2*Nxyzc(&geo)+2, 0, &m->J);
		
	KSPCreate(PETSC_COMM_WORLD,&m->ksp);


	char w[PETSC_MAX_PATH_LEN], filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", m->name, Output_Suffix);
	
	FILE *fp = fopen(filename, "r");
	

	if(fp==NULL){

		sprintf(w, "failed to read %s", filename);
		MyError(w );
	}

   if(GetRank()==0){ fscanf(fp, "%*[^\n]\n", NULL); } // "mode = ["



	ReadVectorC(fp, 2*Nxyzc(&geo)+2, m->vpsi);

   if(GetRank()==0){

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	fgets(w, PETSC_MAX_PATH_LEN, fp);
		
	// Make sure this part is consistent with Mode::Write()
	// make sure to Bcast these at the end

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "ifix=%i\n", &m->ifix); 

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "b=[%i %i %i];", &m->b[0][0], &m->b[1][0], &m->b[2][0]);
	for(int i=0; i<3; i++) m->b[i][1] = 0;

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "BCPeriod=%i;", &m->BCPeriod); 	 

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "D=%lf;", Dout); 	 

	fgets(w, PETSC_MAX_PATH_LEN, fp); // k=[ line
	for(int i=0; i<3; i++){ 
		fgets(w, PETSC_MAX_PATH_LEN, fp);
		sscanf(w, "%lf", &m->k[i]);
	}
   }
   
	fclose(fp);

	MPI_Bcast(m->k, 3, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&m->ifix, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	for(int i=0; i<3; i++) 
	   MPI_Bcast(m->b[i], 2, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&m->BCPeriod, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(Dout, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
		
	double val1, val2;
	GetLast2(m->vpsi, &val1, &val2);
	m->lasing = val2 >= 0.0;

}

void DestroyMode(Mode *m){

	DestroyVec(&m->vpsi);
	DestroyMat(&m->J);
	
	if(!m->ksp){
		KSPDestroy(&m->ksp);
		m->ksp = 0;
	}

}


void Fix(Mode *m, Geometry& geo){

	int N;
	VecGetSize(m->vpsi, &N);
	int Nxyzc = (N-2)/2;

	double max;
	VecMax(m->vpsi,&m->ifix, &max);
	m->ifix = m->ifix % Nxyzc;

	double psifix_real = GetValue(m->vpsi, m->ifix % Nxyzc),
	       psifix_imag = GetValue(m->vpsi, m->ifix % Nxyzc + Nxyzc);

	dcomp factor = OptionsDouble("-norm") / (psifix_real+ ComplexI * psifix_imag);

	VecCopy(m->vpsi, geo.vscratch[0]);
	geo.TimesI(m->vpsi, geo.vscratch[1]);
	
	Complexfun psi;
	CreateComplexfun(&psi, geo.vscratch[0], geo.vscratch[1]);
	Vecfun psiket;
	CreateVecfun(&psiket,m->vpsi);

	for(int i=psiket.ns; i<psiket.ne; i++){

		dcomp val = valc(&psi, i) * factor;
		setr(&psiket, i, ir(&geo, i)? val.imag() : val.real() ); 
	}

	DestroyVecfun(&psiket);

	DestroyComplexfun(&psi);
}



void Write(Mode *m, const Geometry& geo){


	Output(m->vpsi, m->name, "psi");


   if(GetRank()==0){
   
   	

	char filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", m->name, Output_Suffix );

	FILE *fp = fopen(filename, "a");
	fprintf(fp, "ifix=%i;\n", m->ifix);

	fprintf(fp, "b=[%i %i %i];\n", m->b[0][0], m->b[1][0], m->b[2][0]);
	fprintf(fp, "BCPeriod=%i;\n", m->BCPeriod);

	fprintf(fp, "D=%1.15g;\n", geo.D);
	fprintf(fp, "k=[\n%1.15g\n%1.15g\n%1.15g\n];\n", m->k[0], m->k[1], m->k[2]);	
	// "read" constructor for Mode depends on this
	
	fclose(fp);
   }
	MPI_Barrier(PETSC_COMM_WORLD);
	
}


void AddPlaceholders(Mat J, Geometry &geo){

        int ns, ne,N;
        MatGetOwnershipRange(J, &ns, &ne);
        MatGetSize(J, &N, NULL);
        int Nh = N/(Nxyzcr(&geo)+2);

	for(int i=ns; i<ne; i++){ 
	
		if(geo.Last2(i)) continue;		

		for(int jh = 0; jh<Nh; jh++){
			int offset =  jh*(Nxyzcr(&geo)+2);
			
			for(int j=0; j<2; j++) // columns
				MatSetValue(J, i, offset+Nxyzcr(&geo)+j, 0.0, ADD_VALUES);
			
			for(int jr=0; jr<2;jr++) for(int jc=0; jc<geo.Nc; jc++) // tensor
				MatSetValue(J, i, offset+Nxyzc(&geo)*jr + Nxyz(&geo)*jc + i % NJ(&geo) % Nxyz(&geo), 0.0, ADD_VALUES);
				
		}		

	}

	char solver[PETSC_MAX_PATH_LEN];
	OptionsGetString("-pc_factor_mat_solver_package", solver);
	if( 0==strcmp(solver, "pastix" ) ){
		PetscPrintf(PETSC_COMM_WORLD, "pastix detected, symmetrizing nonzero pattern...\n");
	
		MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
		for(int i=ns; i<ne; i++) if( geo.Last2(i) ){
			for(int j = 0; j < N; j++){
				if( geo.Last2(j) ) continue;
				MatSetValue(J, i, j, 0.0, ADD_VALUES);

			}
		}
	}

}

void AddRowDerivatives(Mat J, Geometry& geo, int ifix, int ih){

	int offset =  ih*(Nxyzcr(&geo)+2);

	// last two rows; just put the row derivatives here!
	if(LastProcess()){
		MatSetValue(J, Nxyzcr(&geo)+offset, 
		ifix+offset, 1.0, ADD_VALUES);
	// putting these here takes care of last two residual elements
	// up to scaling of the first one (done below)
		MatSetValue(J, Nxyzcr(&geo) + 1+offset, 
		Nxyzc(&geo)+ ifix+offset, 1.0, ADD_VALUES);
	}	
	

}

void AllocateJacobian(Mat J, Geometry& geo){


	int N, ns, ne;
	MatGetSize(J, &N, NULL);
	MatGetOwnershipRange(J, &ns, &ne);
	int Nh = N / (2*Nxyzc(&geo)+2);
	
	int nnz = 26+Nh*(2+2*geo.Nc);
	// in parenthesis: 2 column derivatives plus 2Nc tensor derivatives

	if(GetSize() > 1) MatMPIAIJSetPreallocation(J, nnz, NULL, nnz, NULL);
	else MatSeqAIJSetPreallocation(J, nnz, NULL);
	// same number for diagonal and off diagonal. Assume N >> size.

}

void Setup(Mode *m, Geometry& geo){


	AllocateJacobian(m->J, geo);
	
    MoperatorGeneralBlochFill(&geo, m->J, m->b, m->BCPeriod, m->k);


	AddPlaceholders(m->J, geo);
	AddRowDerivatives(m->J, geo, m->ifix);
	AssembleMat(m->J);
	MatSetOption(m->J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(m->J); 
	
	KSPSetFromOptions(m->ksp);
	KSPSetOperators(m->ksp, m->J, m->J, SAME_PRECONDITIONER);
		// this will only be called the first time KSPSolve is called. SAME_PRECONDITIONER makes this
		// line analogous to "static". If call it here, the LU factorization
		// will give a better preconditioner than if you call it with just
		// the Moperator, although if you want to share KSPs between different modes
		// with the same BCs you don't want to do this.	

}


double getc(Mode *m){

	if( !m->lasing) return 0.0;
	else return GetFromLast(m->vpsi, 1);

}

dcomp getw(Mode *m){

	if( m->lasing) return GetFromLast(m->vpsi, 0);
	else return GetFromLast(m->vpsi, 0) + ComplexI * GetFromLast(m->vpsi, 1);
}

dcomp gamma_w(Mode *m, Geometry& geo){
	
	return geo.y/( getw(m) -geo.wa + ComplexI*geo.y);
	
}




