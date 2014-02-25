#include <slepc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <new>
#include <cstdio>
#include <complex>
#include <vector>
#include <list>


#include <sys/time.h>
typedef struct timeval tv;

double dt(tv t1, tv t2);

typedef std::complex<double> dcomp;
static const dcomp ComplexI(0.0, 1.0);
static const char Output_Suffix[PETSC_MAX_PATH_LEN] = "_file.m";

int GetRank();
int GetSize();
int LastProcess();
void ScatterRange (Vec x, Vec y, int ix, int iy, int N);
void SetLast2(Vec f, double val1, double val2);
double GetFromLast(Vec v, int ir);
void GetLast2(Vec f, double *val1, double *val2);

void CreateVec(int N, Vec *x);
void AssembleVec(Vec x);
void AssembleMat(Mat M);
void DestroyVec(Vec *x);
void DestroyMat(Mat *A);

int OptionsInt(const char* c);
double OptionsDouble(const char* c);

int OptionsGetInt(const char* c, int* a);
int OptionsGetDouble(const char* c, double* a);
int OptionsGetString(const char* c, char* a);


void OptionsXYZDouble(const char* prefix, double* a);
void OptionsXYZInt(const char* prefix, int* a);



struct Grid{	
	
	Grid(){}
	Grid(int* M, int Mc, int Mr){
		for(int i=0; i<3; i++) N[i] = M[i];
		Nc = Mc;
		Nr = Mr;
	}

	
	int N[3], Nc, Nr;
};


int xyzGrid(Grid *g);
int xyzcGrid(Grid *g);
int xyzcrGrid(Grid *g);


struct Point{

	
	Point(int i, const Grid& H){
		for(int j = 2; j>=0; j--){
			G = H;
			ix[j] = i % G.N[j]; 
			i /= G.N[j];
		}
		ic = i % G.Nc;
		ir = i / G.Nc;	
	}

	Point(int* jx, int jc, int jr, const Grid& H){
		G = H;
		for(int k=0; k<3; k++) ix[k] = jx[k];
		ic = jc;
		ir = jr;
	}

	int ix[3], ir, ic;
	Grid G;
};

int convert(Point *p, int Nc);
int project(Point *p, int Nc);
int projectmedium(Point *p, const Grid& gm, int LowerPML);
int xyz(Point *p);
int xyzc(Point *p);
int xyzcr(Point *p);




#define SCRATCHNUM 5
struct Geometry{

	

	int Npml[3], Nc, LowerPML;
	double h[3];
	Vec vepspml;

	Grid gN, gM;



	Vec vscratch[SCRATCHNUM], vMscratch[SCRATCHNUM], vNhscratch[SCRATCHNUM], vH, veps, vIeps, vf;
	double D, wa, y;
	Geometry();
	~Geometry();	

	void InterpolateVec(Vec vM, Vec vN);
	void CollectVec(Vec vN, Vec vM);
	void Stamp(Vec vN, int ic, int ir, Vec scratchM);
	void TimesI(Vec v, Vec Iv);
	void VecSqMedium(Vec v, Vec vsq, Vec scratchM);

	int Last2(int i);
	void SetJacobian(Mat J, Vec v, int jc, int jr, int jh);

};

void MoperatorGeneralBlochFill(Geometry *geo, Mat A,  int b[3][2], int DimPeriod, double k[3], int ih=0);	

int Nxyz(Geometry *geo);
int Nxyzc(Geometry *geo);
int Nxyzcr(Geometry *geo);
int Mxyz(Geometry *geo);

int NJ(Geometry *geo);
int offset(Geometry *geo, int ih);
int ir(Geometry *geo, int i);

struct Mode{

	

	Vec vpsi;
	char name[PETSC_MAX_PATH_LEN];
	double k[3];
	int ifix, BCPeriod, b[3][2], lasing;
	Mat J;
	KSP ksp; // one ksp per J seems faster

	Mode(Geometry& geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]);
	Mode(char *Name, Geometry& geo, double *Dout); // read constructor
	~Mode();
	void Fix(Geometry& geo);
	void Setup(Geometry& geo);
	void Write(const Geometry& geo);
	double c();
	dcomp w();
	dcomp gamma_w(Geometry& geo);

};
typedef std::list<Mode*> modelist;
#define FORMODES(L, it)   for(modelist::iterator it=L.begin(); it!= L.end(); it++) 
// note no ; at end of macro!

void ComputeGain(Geometry& geo, modelist& L);
// not sure how to define this as a member function of Geometry, since Mode is defined after geometry

void CreateSquareMatrix(int N, int nz, Mat *A);
double GetValue(Vec v, int i);
void ReadVectorC(FILE *fp, int N, Vec v);

PetscErrorCode MyError(const char* message);
double sqr(const double a);
dcomp sqr(dcomp a);
void View(Vec x, PetscViewer viewer);

void Output(Vec A, const char* name, const char* variable_name);


void NewtonSolve(modelist &L, Geometry& geo, Vec v, Vec f, Vec dv);
void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, modelist &Lh, Vec vNh, Mode& m, Geometry& geo, Vec f, Vec dv);
double FormJf(modelist& L, Geometry& geo, Vec v, Vec f);


struct Vecfun{
// always Nxyzcr()+2!	
	
	Vec u;
	int ms, me; // should keep these, since they are used by the val() function
				// every single time.

	 
	double* a;

	
	Vecfun(Vec w){
		u = w;
		VecGetOwnershipRange(u, &ms, &me); // TODO: take into account Nxyzcr+2 here
		int N;
		VecGetSize(w, &N);
		if(me >= N-2) me = N-2;
		
		VecGetArray(u, &a);
	}
	~Vecfun(){ VecRestoreArray(u, &a); }

	int ns() const{ return ms;}
	int ne() const{ return me;}
	double valr(int i) const{ return a[i-ns()];}
	void set(int i, double val){ a[i-ns()] = val;}

};

struct ComplexVecfun: public Vecfun{

	
	Vec v;
	double* b;
	int Nxyzc;
	
	
	ComplexVecfun(Vec w, Vec x): Vecfun(w){
		v = x;
		VecGetArray(v, &b);
		
		int N;
		VecGetSize(v, &N);
		Nxyzc = (N-2)/2;
	}

	void set(int i, dcomp val){
		
		if(i/Nxyzc) val /= ComplexI;
		a[i-ns()] = val.real();
		b[i-ns()] = -val.imag();
	}
	~ComplexVecfun(){ VecRestoreArray(v, &b); }

};

dcomp valc(ComplexVecfun *fun, int i);


dcomp pmlval(int i, int* N, int* Npml, double* h, int LowerPML, int k);
void AddPlaceholders(Mat J, Geometry &geo);
void AllocateJacobian(Mat J, Geometry& geo);
void AddRowDerivatives(Mat J, Geometry& geo, int ifix, int ih=0);



