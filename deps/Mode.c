#include "headers.h"

Mode CreateMode(Geometry geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]){
	Mode m = (Mode) malloc(sizeof(struct Mode_s) );
	CreateVec(2*Nxyzc(geo)+2, &m->vpsi);
	CreateSquareMatrix(2*Nxyzc(geo)+2, 0, &m->J);
	KSPCreate(PETSC_COMM_WORLD,&m->ksp);

	m->lasing = 0;
	m->ifix = ifix_;
	m->BCPeriod = BCPeriod_;
	int i, j;
	for(i=0; i<3; i++) for(j=0; j<2; j++) m->b[i][j] = b_[i][j];
	for(i=0; i<3; i++) m->k[i] = k_[i];

	return m;
}

Mode ModeRead(const char *Name, Geometry geo, double *Dout){
// TODO free after every ModeRead

	Mode m = (Mode) malloc(sizeof(struct Mode_s) );
	sprintf(m->name, "%s", Name );
	CreateVec(2*Nxyzc(geo)+2, &m->vpsi);
	CreateSquareMatrix(2*Nxyzc(geo)+2, 0, &m->J);

	KSPCreate(PETSC_COMM_WORLD,&m->ksp);

	char w[PETSC_MAX_PATH_LEN], filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", m->name, Output_Suffix);
		int i;
	FILE *fp = fopen(filename, "r");

	if(fp==NULL){
		sprintf(w, "failed to read %s", filename);
		MyError(w );
	}

   if(GetRank()==0){ fgets(w, PETSC_MAX_PATH_LEN, fp); } // "mode = ["

	ReadVectorC(fp, 2*Nxyzc(geo)+2, m->vpsi);

   if(GetRank()==0){
	fgets(w, PETSC_MAX_PATH_LEN, fp);
	fgets(w, PETSC_MAX_PATH_LEN, fp);

	// Make sure this part is consistent with Mode::Write()
	// make sure to Bcast these at the end

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "ifix=%i\n", &m->ifix); 

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "b=[%i %i %i];", &m->b[0][0], &m->b[1][0], &m->b[2][0]);

	for(i=0; i<3; i++) m->b[i][1] = 0;

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "BCPeriod=%i;", &m->BCPeriod); 	 

	fgets(w, PETSC_MAX_PATH_LEN, fp);
	sscanf(w, "D=%lf;", Dout); 	 

	fgets(w, PETSC_MAX_PATH_LEN, fp); // k=[ line
	for(i=0; i<3; i++){ 
		fgets(w, PETSC_MAX_PATH_LEN, fp);
		sscanf(w, "%lf", &m->k[i]);
	}
   }
   
	fclose(fp);

	MPI_Bcast(m->k, 3, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&m->ifix, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	for(i=0; i<3; i++) 
	   MPI_Bcast(m->b[i], 2, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&m->BCPeriod, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(Dout, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

	double val1, val2;
	GetLast2(m->vpsi, &val1, &val2);
	m->lasing = val2 >= 0.0;

	return m;
}

void DestroyMode(Mode m){
	DestroyVec(&m->vpsi);
	DestroyMat(&m->J);

	if(!m->ksp){
		KSPDestroy(&m->ksp);
		m->ksp = 0;
	}
}

void CopyPsi(Mode m, double *psiout){
	int ns, ne;
	VecGetOwnershipRange(m->vpsi, &ns, &ne);

	double *psi;
	VecGetArray(m->vpsi, &psi);
	int i;
	for(i=ns; i<ne; i++)
		psiout[i-ns] = psi[i-ns];

	VecRestoreArray(m->vpsi, &psi);
}

int PsiSize(Mode m){
	int N;
	VecGetSize(m->vpsi, &N);
	return N;
}

void Fix(Mode m, Geometry geo, double norm){
	int N;
	VecGetSize(m->vpsi, &N);
	int Nxyzc = (N-2)/2, i;

	double max;
	VecMax(m->vpsi,&m->ifix, &max);
	m->ifix = m->ifix % Nxyzc;

	double psifix_real = GetValue(m->vpsi, m->ifix % Nxyzc),
	       psifix_imag = GetValue(m->vpsi, m->ifix % Nxyzc + Nxyzc);

	dcomp factor = norm / (psifix_real+ ComplexI * psifix_imag);

	VecCopy(m->vpsi, geo->vscratch[0]);
	TimesI(geo, m->vpsi, geo->vscratch[1]);

	Complexfun psi;
	CreateComplexfun(&psi, geo->vscratch[0], geo->vscratch[1]);
	Vecfun psiket;
	CreateVecfun(&psiket,m->vpsi);

	for(i=psiket.ns; i<psiket.ne; i++){
		dcomp val = valc(&psi, i) * factor;
		setr(&psiket, i, ir(geo, i)? cimag(val) : creal(val) ); 
	}

	DestroyVecfun(&psiket);

	DestroyComplexfun(&psi);
}

void Write(Mode m, const Geometry geo){
	Output(m->vpsi, m->name, "psi");

   if(GetRank()==0){
   
   

	char filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", m->name, Output_Suffix );

	FILE *fp = fopen(filename, "a");
	fprintf(fp, "ifix=%i;\n", m->ifix);

	fprintf(fp, "b=[%i %i %i];\n", m->b[0][0], m->b[1][0], m->b[2][0]);
	fprintf(fp, "BCPeriod=%i;\n", m->BCPeriod);

	fprintf(fp, "D=%1.15g;\n", geo->D);
	fprintf(fp, "k=[\n%1.15g\n%1.15g\n%1.15g\n];\n", m->k[0], m->k[1], m->k[2]);
	// "read" constructor for Mode depends on this

	// additional lines for plotting only, not read in ReadMode
	const int *N = &(geo->gN.N[0]);
	const double *h = &(geo->h[0]);
	fprintf(fp, "N = [%i, %i, %i];\n", N[0],  N[1], N[2]);
	fprintf(fp, "h = [%1.8g, %1.8g, %1.8g];\n", h[0],  h[1], h[2]);
	fprintf(fp, "LowerPML = %i;\n", geo->LowerPML);
	fprintf(fp, "Nc = %i;\n", geo->Nc);

	fclose(fp);
   }
	MPI_Barrier(PETSC_COMM_WORLD);
}

void AddRowDerivatives(Mat J, Geometry geo, int ifix, int ih){
	int offset =  ih*(Nxyzcr(geo)+2);

	// last two rows; just put the row derivatives here!
	if(LastProcess()){
		MatSetValue(J, Nxyzcr(geo)+offset, 
		ifix+offset, 1.0, ADD_VALUES);
	// putting these here takes care of last two residual elements
	// up to scaling of the first one (done below)
		MatSetValue(J, Nxyzcr(geo) + 1+offset, 
		Nxyzc(geo)+ ifix+offset, 1.0, ADD_VALUES);
	}
}

void AllocateJacobian(Mat J, Geometry geo){
	int N, ns, ne;
	MatGetSize(J, &N, NULL);
	MatGetOwnershipRange(J, &ns, &ne);
	int Nh = N / NJ(geo);

	int nnz = 26+Nh*(2+2*geo->Nc);
	// in parenthesis: 2 column derivatives plus 2Nc tensor derivatives

	if(GetSize() > 1) MatMPIAIJSetPreallocation(J, nnz, NULL, nnz, NULL);
	else MatSeqAIJSetPreallocation(J, nnz, NULL);
	// same number for diagonal and off diagonal. Assume N >> size.

	//================== Adding placeholders ================

        int i, jh, j, jr, jc;
        MatGetOwnershipRange(J, &ns, &ne);
        MatGetSize(J, &N, NULL);

	for(i=ns; i<ne; i++){ 

		if(!Last2(geo, i)) for(jh = 0; jh<Nh; jh++){
			int offset =  jh*(Nxyzcr(geo)+2);

			for(j=0; j<2; j++) // columns
				MatSetValue(J, i, offset+Nxyzcr(geo)+j, 0.0, ADD_VALUES);

			for(jr=0; jr<2;jr++) for(jc=0; jc<geo->Nc; jc++) // tensor
				MatSetValue(J, i, offset+Nxyzc(geo)*jr + Nxyz(geo)*jc + i % NJ(geo) % Nxyz(geo), 0.0, ADD_VALUES);
		}
	}

	/*    // this works, just disabling for interface reasons
	char solver[PETSC_MAX_PATH_LEN];
	OptionsGetString("-pc_factor_mat_solver_package", solver);
	if( 0==strcmp(solver, "pastix" ) ){
		PetscPrintf(PETSC_COMM_WORLD, "pastix detected, symmetrizing nonzero pattern...\n");

		MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
		for(i=ns; i<ne; i++) if( Last2(geo, i) ){
			for(j = 0; j < N; j++){
				if( Last2(geo, j) ) continue;
				MatSetValue(J, i, j, 0.0, ADD_VALUES);

			}
		}
	}
	*/
}

void Setup(Mode m, Geometry geo){
	AllocateJacobian(m->J, geo);

    MoperatorGeneralBlochFill(geo, m->J, m->b, m->BCPeriod, m->k, 0);
	AddRowDerivatives(m->J, geo, m->ifix, 0);
	AssembleMat(m->J);
	MatSetOption(m->J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(m->J); 

	KSPSetFromOptions(m->ksp);
	PC pc;
	KSPGetPC(m->ksp,&pc);
 	PCSetType(pc,PCLU);
  	PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
	// don't forget to change this in Bundle too

	KSPSetOperators(m->ksp, m->J, m->J, SAME_PRECONDITIONER);
		// this will only be called the first time KSPSolve is called. SAME_PRECONDITIONER makes this
		// line analogous to "static". If call it here, the LU factorization
		// will give a better preconditioner than if you call it with just
		// the Moperator, although if you want to share KSPs between different modes
		// with the same BCs you don't want to do this.
}

double get_c(Mode m){
	if( !m->lasing) return 0.0;
	else return GetFromLast(m->vpsi, 1);
}

dcomp get_w(Mode m){
	if( m->lasing) return GetFromLast(m->vpsi, 0);
	else return GetFromLast(m->vpsi, 0) + ComplexI * GetFromLast(m->vpsi, 1);
}

dcomp gamma_w(Mode m, Geometry geo){
	return geo->y/( get_w(m) -geo->wa + ComplexI*geo->y);
}

void addArrayMode(Mode **ma, int old_size, Mode m){
	if(old_size == 0){
		*ma = (Mode *) malloc( sizeof(struct Mode_s) );
	}else{
		*ma = (Mode *) realloc( *ma, (old_size+1)*sizeof(struct Mode_s) );
	}

	(*ma)[old_size] = m;
}

Mode GetMode(Mode *ms, int n){ return ms[n] ; }

