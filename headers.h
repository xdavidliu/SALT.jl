#include <slepc.h>
#include <cmath>
#include <cstdio>
#include <complex>
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



typedef struct Grid_s{	
	
	int N[3], Nc, Nr;
} Grid;


void CreateGrid(Grid *g, int* M, int Mc, int Mr);


int xyzGrid(Grid *g);
int xyzcGrid(Grid *g);
int xyzcrGrid(Grid *g);


typedef struct Point_s{


	int ix[3], ir, ic;
	Grid G;
} Point;

void CreatePoint_i(Point *p, int i, Grid *H);

int convert(Point *p, int Nc);
int project(Point *p, int Nc);
int projectmedium(Point *p, Grid *gm, int LowerPML);
int xyz(Point *p);
int xyzc(Point *p);
int xyzcr(Point *p);




#define SCRATCHNUM 5
typedef struct Geometry_s{

	

	int Npml[3], Nc, LowerPML;
	double h[3];
	Vec vepspml;

	Grid gN, gM;



	Vec vscratch[SCRATCHNUM], vMscratch[SCRATCHNUM], vNhscratch[SCRATCHNUM], vH, veps, vIeps, vf;
	double D, wa, y;

} Geometry;

void CreateGeometry(Geometry *geo);
void DestroyGeometry(Geometry *geo);	
void InterpolateVec(Geometry *geo, Vec vM, Vec vN);
void CollectVec(Geometry *geo, Vec vN, Vec vM);
void Stamp(Geometry *geo, Vec vN, int ic, int ir, Vec scratchM);
void TimesI(Geometry *geo, Vec v, Vec Iv);
void VecSqMedium(Geometry *geo, Vec v, Vec vsq, Vec scratchM);
int Last2(Geometry *geo, int i);
void SetJacobian(Geometry *geo, Mat J, Vec v, int jc, int jr, int jh);

void MoperatorGeneralBlochFill(Geometry *geo, Mat A,  int b[3][2], int DimPeriod, double k[3], int ih=0);	

int Nxyz(Geometry *geo);
int Nxyzc(Geometry *geo);
int Nxyzcr(Geometry *geo);
int Mxyz(Geometry *geo);

int NJ(Geometry *geo);
int offset(Geometry *geo, int ih);
int ir(Geometry *geo, int i);

typedef struct Mode_s{

	Vec vpsi;
	char name[PETSC_MAX_PATH_LEN];
	double k[3];
	int ifix, BCPeriod, b[3][2], lasing;
	Mat J;
	KSP ksp; // one ksp per J seems faster

} Mode;

void CreateMode(Mode *m, Geometry *geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]);
void ModeRead(Mode *m, char *Name, Geometry *geo, double *Dout);
void DestroyMode(Mode *m);
void Setup(Mode *m, Geometry *geo);

void Write(Mode *m, const Geometry *geo);
double getc(Mode *m);
dcomp getw(Mode *m);
dcomp gamma_w(Mode *m, Geometry *geo);
void Fix(Mode *m, Geometry *geo);

typedef std::list<Mode*> modelist;


typedef struct ModeArray_s{
	int size;
	Mode **L;
} ModeArray;

void CreateModeArray(ModeArray *ma, Mode *m);
void DestroyModeArray(ModeArray *ma);
void AddArrayMode(ModeArray *ma, Mode *m);
void RemoveArrayMode(ModeArray *ma, int n);

void CreateFromList(ModeArray *ma, modelist& L); // temp

#define FORMODES(L, it)   for(modelist::iterator it=L.begin(); it!= L.end(); it++) 
// note no ; at end of macro!

void ComputeGain(Geometry *geo, modelist& L);
// not sure how to define this as a member function of Geometry, since Mode is defined after geometry

void CreateSquareMatrix(int N, int nz, Mat *A);
double GetValue(Vec v, int i);
void ReadVectorC(FILE *fp, int N, Vec v);

PetscErrorCode MyError(const char* message);
double sqr(const double a);
dcomp sqr(dcomp a);
void View(Vec x, PetscViewer viewer);

void Output(Vec A, const char* name, const char* variable_name);


void NewtonSolve(modelist &L, Geometry *geo, Vec v, Vec f, Vec dv);
void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, modelist &Lh, Vec vNh, Mode *m, Geometry *geo, Vec f, Vec dv);
double FormJf(modelist& L, Geometry *geo, Vec v, Vec f);


typedef struct Vecfun_s{
// always Nxyzcr()+2!	
	
	Vec u;
	int ns, ne; // should keep these, since they are used by the val() function
				// every single time.

	 
	double* a;


} Vecfun;

void CreateVecfun(Vecfun *fun, Vec w);

double valr(Vecfun *fun, int i);
void setr(Vecfun *fun, int i, double val);
void DestroyVecfun(Vecfun *fun);


typedef struct Complexfun_s{

	Vec u, v;
	double *a, *b;
	int Nxyzc, ns, ne;

}Complexfun;

dcomp valc(Complexfun *fun, int i);
void setc(Complexfun *fun,int i, dcomp val);

void CreateComplexfun(Complexfun *fun, Vec w, Vec x);
void DestroyComplexfun(Complexfun *fun);

dcomp pmlval(int i, int* N, int* Npml, double* h, int LowerPML, int k);
void AddPlaceholders(Mat J, Geometry *geo);
void AllocateJacobian(Mat J, Geometry *geo);
void AddRowDerivatives(Mat J, Geometry *geo, int ifix, int ih=0);



