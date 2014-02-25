#include "headers.h"

double dt(tv t1, tv t2){

	double val = (t2.tv_sec - t1.tv_sec)*1000
			+ (t2.tv_usec - t1.tv_usec)/1000.0;
	return val/1000;
}


PetscErrorCode MyError(const char* message){

	SETERRQ(PETSC_COMM_WORLD, 1, message);
}

double sqr(double a){ return a*a;}
dcomp sqr(dcomp a){ return a*a;}

void AssembleVec(Vec x){ VecAssemblyBegin(x); VecAssemblyEnd(x); }
void AssembleMat(Mat M){ MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);}
void DestroyVec(Vec *x){
	if(!x) VecDestroy(x);
	x = 0;	
}

void DestroyMat(Mat *A){
	if(!A) MatDestroy(A); 
	A = 0;
}

void View(Mat A, PetscViewer viewer){ MatView(A, viewer); }
void View(Vec x, PetscViewer viewer){ VecView(x, viewer); }

void OptionsXYZDouble(const char* prefix, double* a){

	char option[PETSC_MAX_PATH_LEN];
	const char x[3] = {'x', 'y', 'z'};
	for(int i=0; i<3; i++){
		sprintf(option, "%s%c", prefix, x[i]);
		OptionsGetDouble(option, &a[i]);
	}

}


void OptionsXYZInt(const char* prefix, int* a){

	char option[PETSC_MAX_PATH_LEN];
	const char x[3] = {'x', 'y', 'z'};
	for(int i=0; i<3; i++){
		sprintf(option, "%s%c", prefix, x[i]);
		OptionsGetInt(option, &a[i]);
	}

}

int OptionsGetString(const char* c, char* a){ 
	PetscBool flg;
	PetscOptionsGetString(PETSC_NULL,c, a, PETSC_MAX_PATH_LEN, &flg); 
	return flg;
}

int OptionsGetInt(const char* c, int* a){
	PetscBool flg;
	PetscOptionsGetInt(PETSC_NULL,c,a,&flg);
	return flg;
}

int OptionsInt(const char* c){
	int out;
	OptionsGetInt(c, &out);
	return out;
}

std::string OptionsString(const char *c){

	char a[PETSC_MAX_PATH_LEN];
	OptionsGetString(c, a);
	return std::string(a);

}

double OptionsDouble(const char* c){
	double out;
	OptionsGetDouble(c, &out);
	return out;
}

int OptionsGetDouble(const char* c, double* a){
	PetscBool flg;
	PetscOptionsGetReal(PETSC_NULL,c,a,&flg); 
	return flg;
}


void ScatterRange (Vec x, Vec y, int ix, int iy, int N){

	int ns, ne;
	VecGetOwnershipRange(x, &ns, &ne);

	const double *vals;
	VecGetArrayRead(x, &vals);

	int i=ns;
	if( ns < ix) i = ix; // don't modify ns because need it to access vals
	if( ne > ix+N) ne = ix+N;

	for(; i<ne; i++)
		VecSetValue(y, i+iy-ix, vals[i-ns], INSERT_VALUES);

	
	VecRestoreArrayRead(x, &vals);
	AssembleVec(y);

	
}


void CreateVec(int N, Vec *x){
	VecCreate(PETSC_COMM_WORLD, x);
	VecSetSizes(*x, PETSC_DECIDE, N);
	VecSetFromOptions(*x);
}

void CreateSquareMatrix(int N, int nz, Mat *A){
	
	if(GetSize() == 1) MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, nz, NULL, A);
	else 	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
		N, N, nz, NULL, nz, NULL, A);
}

double GetValue(Vec v, int i){

	const int *ranges; int owner=0;
	VecGetOwnershipRanges(v, &ranges);
	for(int j=0; j<GetSize(); j++){
		if( ranges[j] <= i && i < ranges[j+1]){
			owner = j;
			break;
		}
	}

	double val;
	if( GetRank() == owner)
		VecGetValues(v, 1, &i, &val);

	MPI_Bcast(&val, 1, MPI_DOUBLE, owner, PETSC_COMM_WORLD);
	return val;
}

void ReadVector(std::ifstream& is, int N, Vec v){

   double val;

   if(GetRank() ==0){
	for(int i=0; i< N; i++){
		is >> val;
		VecSetValue(v, i, val, INSERT_VALUES);
	}
   }
   AssembleVec(v);
}



void SetLast2(Vec f, double val1, double val2){

	int N;
	VecGetSize(f, &N);

	if(LastProcess()){
		VecSetValue(f, N-2, val1, INSERT_VALUES);
		VecSetValue(f, N-1, val2, INSERT_VALUES);
	}
	AssembleVec(f);

}

double GetFromLast(Vec v, int ir){
// MUCH faster than using single element scatters!

	int N, row;
	VecGetSize(v, &N);
	row = N - 2+ir;
	
	double val;
	if( LastProcess() )
		VecGetValues(v, 1, &row, &val);
		
	MPI_Bcast(&val, 1, MPI_DOUBLE, GetSize()-1, PETSC_COMM_WORLD);
	return val;
}

void GetLast2(Vec f, double *val1, double *val2){

	int N;
	VecGetSize(f, &N);
	if(val1!=NULL) *val1 = GetFromLast(f, 0);
	if(val2!=NULL) *val2 = GetFromLast(f, 1);	
}

int GetRank(){
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	return rank;
}

int GetSize(){
	int size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	return size;
}

int LastProcess(){ 
	return GetRank() == GetSize()-1;
}

void Output(Vec A, const char* name, const char* variable_name){

	char filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", name, Output_Suffix.c_str());

	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
	PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject) A, strncmp(variable_name, "", PETSC_MAX_PATH_LEN) ?variable_name : name);
	View(A, viewer);
	PetscViewerDestroy(&viewer);
}