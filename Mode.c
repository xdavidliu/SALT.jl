#include "headers.h"




Mode::Mode(Geometry& geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]){

	CreateVec(2*geo.Nxyzc()+2, &vpsi);
	CreateSquareMatrix(2*geo.Nxyzc()+2, 0, &J);
	KSPCreate(PETSC_COMM_WORLD,&ksp);

	lasing = 0;
	ifix = ifix_;
	BCPeriod = BCPeriod_;
	for(int i=0; i<3; i++) for(int j=0; j<2; j++) b[i][j] = b_[i][j];
	for(int i=0; i<3; i++) k[i] = k_[i];
}



Mode::Mode(std::string Name, Geometry& geo, double *Dout){ // read constructor

	name = Name;
	CreateVec(2*geo.Nxyzc()+2, &vpsi);
	CreateSquareMatrix(2*geo.Nxyzc()+2, 0, &J);
		
	KSPCreate(PETSC_COMM_WORLD,&ksp);




	std::string filename = name + Output_Suffix;
	std::ifstream is(filename.c_str());
	std::string s;

	if(!is){
		s = "failed to read ";
		s += filename;
		MyError(s.c_str() );
	}

   if(GetRank()==0){ std::getline(is, s); } // "mode = ["



	ReadVector(is, 2*geo.Nxyzc()+2, vpsi);

   if(GetRank()==0){
	std::getline(is, s);
	std::getline(is, s); // for some reason you need two getlines to get the "];" character here

	// Make sure this part is consistent with Mode::Write()
	// make sure to Bcast these at the end

	std::getline(is, s);
	std::sscanf(s.c_str(), "ifix=%i;", &ifix); 

	std::getline(is, s);
	std::sscanf(s.c_str(), "b=[%i %i %i];", 
		&b[0][0], &b[1][0], &b[2][0]);
	for(int i=0; i<3; i++) b[i][1] = 0;

	std::getline(is, s);
	std::sscanf(s.c_str(), "BCPeriod=%i;", &BCPeriod); 	 

	std::getline(is, s);
	std::sscanf(s.c_str(), "D=%lf;", Dout); 	 

	std::getline(is, s); // k=[ line
	for(int i=0; i<3; i++){ 
		std::getline(is, s);
		std::sscanf(s.c_str(), "%lf", &k[i]);
	}
   }

	MPI_Bcast(k, 3, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&ifix, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	for(int i=0; i<3; i++) 
	   MPI_Bcast(b[i], 2, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&BCPeriod, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(Dout, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
		
	double val1, val2;
	GetLast2(vpsi, &val1, &val2);
	lasing = val2 >= 0.0;

}

Mode::~Mode(){

	Destroy(&vpsi);
	Destroy(&J);
	
	if(!ksp){
		KSPDestroy(&ksp);
		ksp = 0;
	}

}


void Mode::Fix(Geometry& geo){

	int N;
	VecGetSize(vpsi, &N);
	int Nxyzc = (N-2)/2;

	double max;
	VecMax(vpsi,&ifix, &max);
	ifix = ifix % Nxyzc;

	double psifix_real = GetValue(vpsi, ifix % Nxyzc),
	       psifix_imag = GetValue(vpsi, ifix % Nxyzc + Nxyzc);

	dcomp factor = OptionsDouble("-norm") / (psifix_real+ ComplexI * psifix_imag);

	VecCopy(vpsi, geo.vscratch[0]);
	geo.TimesI(vpsi, geo.vscratch[1]);
	
	{Vecfun psiket(vpsi); ComplexVecfun psi(geo.vscratch[0], geo.vscratch[1]);

	for(int i=psiket.ns(); i<psiket.ne(); i++){

		dcomp val = psi.val(i) * factor;
		psiket.set(i, geo.ir(i)? val.imag() : val.real() ); 
	}}

}



void Mode::Write(const Geometry& geo){


	Output(vpsi, name.c_str(), "psi");


   if(GetRank()==0){
	std::string filename = name + Output_Suffix;
	std::ofstream f(filename.c_str(), std::ios_base::app);

	f << "ifix=" << ifix << ";\n";
	f << "b=[" << b[0][0] << " " << b[1][0] << " "
	  << b[2][0] << "];\n";
	f << "BCPeriod=" << BCPeriod << ";\n";
	f << std::setprecision(15) << "D=" << geo.D << ";\n";
	f << "k=[\n" << k[0] << "\n" << k[1] << "\n" << k[2] << "\n];\n";
	// "read" constructor for Mode depends on this
   }
	MPI_Barrier(PETSC_COMM_WORLD);
	
}


void AddPlaceholders(Mat J, Geometry &geo){

        int ns, ne,N;
        MatGetOwnershipRange(J, &ns, &ne);
        MatGetSize(J, &N, NULL);
        int Nh = N/(geo.Nxyzcr()+2);

	for(int i=ns; i<ne; i++){ 
	
		if(geo.Last2(i)) continue;		

		for(int jh = 0; jh<Nh; jh++){
			int offset =  jh*(geo.Nxyzcr()+2);
			
			for(int j=0; j<2; j++) // columns
				MatSetValue(J, i, offset+geo.Nxyzcr()+j, 0.0, ADD_VALUES);
			
			for(int jr=0; jr<2;jr++) for(int jc=0; jc<geo.Nc; jc++) // tensor
				MatSetValue(J, i, offset+geo.Nxyzc()*jr + geo.Nxyz()*jc + i % geo.NJ() % geo.Nxyz(), 0.0, ADD_VALUES);
				
		}		

	}

	if( 0==strcmp(OptionsString("-pc_factor_mat_solver_package").c_str(), "pastix" ) ){
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

	int offset =  ih*(geo.Nxyzcr()+2);

	// last two rows; just put the row derivatives here!
	if(LastProcess()){
		MatSetValue(J, geo.Nxyzcr()+offset, 
		ifix+offset, 1.0, ADD_VALUES);
	// putting these here takes care of last two residual elements
	// up to scaling of the first one (done below)
		MatSetValue(J, geo.Nxyzcr() + 1+offset, 
		geo.Nxyzc()+ ifix+offset, 1.0, ADD_VALUES);
	}	
	

}

void AllocateJacobian(Mat J, Geometry& geo){


	int N, ns, ne;
	MatGetSize(J, &N, NULL);
	MatGetOwnershipRange(J, &ns, &ne);
	int Nh = N / (2*geo.Nxyzc()+2);
	
	int nnz = 26+Nh*(2+2*geo.Nc);
	// in parenthesis: 2 column derivatives plus 2Nc tensor derivatives

	if(GetSize() > 1) MatMPIAIJSetPreallocation(J, nnz, NULL, nnz, NULL);
	else MatSeqAIJSetPreallocation(J, nnz, NULL);
	// same number for diagonal and off diagonal. Assume N >> size.

}

void Mode::Setup(Geometry& geo){


	AllocateJacobian(J, geo);
	
        geo.MoperatorGeneralBlochFill(J, b, BCPeriod, k);


	AddPlaceholders(J, geo);
	AddRowDerivatives(J, geo, ifix);
	Assemble(J);
	MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(J); 
	
	KSPSetFromOptions(ksp);
	KSPSetOperators(ksp, J, J, SAME_PRECONDITIONER);
		// this will only be called the first time KSPSolve is called. SAME_PRECONDITIONER makes this
		// line analogous to "static". If call it here, the LU factorization
		// will give a better preconditioner than if you call it with just
		// the Moperator, although if you want to share KSPs between different modes
		// with the same BCs you don't want to do this.	

}


double Mode::c(){

	if( !lasing) return 0.0;
	else return GetFromLast(vpsi, 1);

}

dcomp Mode::w(){

	if( lasing) return GetFromLast(vpsi, 0);
	else return GetFromLast(vpsi, 0) + ComplexI * GetFromLast(vpsi, 1);
}

dcomp Mode::gamma_w(Geometry& geo){
	
	return geo.y/( w() -geo.wa + ComplexI*geo.y);
	
}




