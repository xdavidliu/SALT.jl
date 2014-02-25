#include <slepc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <new>
#include <cstdio>
#include <complex>
#include <string>
#include <vector>
#include <list>


#include <sys/time.h>
typedef struct timeval tv;

double dt(tv t1, tv t2);

typedef std::complex<double> dcomp;
static const dcomp ComplexI(0.0, 1.0);
const std::string Output_Suffix = "_file.m";

int GetRank();
int GetSize();
int LastProcess();
void ScatterRange (Vec x, Vec y, int ix, int iy, int N);
void SetLast2(Vec f, double val1, double val2);
double GetFromLast(Vec v, int ir);
void GetLast2(Vec f, double *val1, double *val2);

void CreateVec(int N, Vec *x);
void Assemble(Vec x);
void Assemble(Mat M);
void Destroy(Vec *x);
void Destroy(Mat *A);

int OptionsInt(const char* c);
double OptionsDouble(const char* c);
std::string OptionsString(const char *c);

int OptionsGet(const char* c, int* a);
int OptionsGet(const char* c, double* a);
int OptionsGet(const char* c, char* a);


template<class T> void OptionsXYZ(const char* prefix, T* a){

	char option[PETSC_MAX_PATH_LEN];
	const char x[3] = {'x', 'y', 'z'};
	for(int i=0; i<3; i++){
		sprintf(option, "%s%c", prefix, x[i]);
		OptionsGet(option, &a[i]);
	}

}



class Grid{	
	public:
	Grid(){}
	Grid(int* M, int Mc, int Mr){
		for(int i=0; i<3; i++) N[i] = M[i];
		Nc = Mc;
		Nr = Mr;
	}
	int x(int ic) const{return N[ic];}
	int c() const{return Nc;}
	int r() const{return Nr;}
	int xyz() const{return N[0]*N[1]*N[2];}
	int xyzc() const{return xyz()*Nc;}
	int xyzcr() const{return xyzc()*Nr;}

	void setc(int Mc){ Nc = Mc;}

	private:
	int N[3], Nc, Nr;
};




class Point{

	public:
	Point(int i, const Grid& H){
		for(int j = 2; j>=0; j--){
			G = H;
			ix[j] = i % G.x(j); 
			i /= G.x(j);
		}
		ic = i % G.c();
		ir = i / G.c();	
	}

	Point(int* jx, int jc, int jr, const Grid& H){
		G = H;
		for(int k=0; k<3; k++) ix[k] = jx[k];
		ic = jc;
		ir = jr;
	}


	int x(int jc) const{return ix[jc];}
	int c() const{return ic;}
	int r() const{return ir;}
	int xyz() const{return x(0)*G.x(2)*G.x(1) + x(1)*G.x(2) + x(2);}
	int xyzc() const{return ic*G.xyz() + xyz();}
	int xyzcr() const{return ir*G.xyzc() + ic*G.xyz() + xyz();}

	void setc(int jc){ ic = jc;}
	void setr(int jr){ ir = jr;}

	int convert(int Nc){

		if(Nc==G.c() ) return ic;
		else if(Nc == 3){
			if(G.c()==1 && ic == 0) return 2;  // TM to vector
			else if(G.c()==2 && ic < 2) return ic;  // TE to vector
			else return -1;
		}else if(G.c() == 3){
			if(Nc==1 && ic == 2) return 0;  // vector to TM
			else if(Nc==2 && ic < 2) return ic;   // vector to TE 
			else return -1;
		}else return -1;

	}


	int project(int Nc){
		ic = convert(Nc);
		G.setc(Nc);
		return ic;
	}


	int projectmedium(const Grid& gm, int LowerPML){

		int medium =1;
		for(int j=0; j<3; j++){ // position component

			double d = x(j) - LowerPML*floor( (G.x(j)-gm.x(j))/2.0 ) + ( c()!=j)*0.5;
			ix[j] = ceil(d-0.5);
			if(x(j)<0 || x(j)>= gm.x(j) ) medium = 0;
		}
		G = gm;
		return medium;
	}

	private:
	int ix[3], ir, ic;
	Grid G;
};






#define SCRATCHNUM 5
class Geometry{

	public:

	int Npml[3], Nc, LowerPML;
	double h[3];
	Vec vepspml;

	Grid gN, gM;

	int Nxyz(){ return gN.xyz(); }
	int Nxyzc(){ return gN.xyzc(); }
	int Nxyzcr(){ return gN.xyzcr();}
	int Mxyz(){ return gM.xyz(); }
	int NJ(){ return Nxyzcr() + 2;}
	int offset(int ih){ return ih*NJ(); }
	int ir(int i){ return i%NJ() / Nxyzc(); }

	void MoperatorGeneralBlochFill(Mat A,  int b[3][2], int DimPeriod, double k[3], int ih=0);	

	Vec vscratch[SCRATCHNUM], vMscratch[SCRATCHNUM], vNhscratch[SCRATCHNUM], vH, veps, vIeps, vf;
	double D, wa, y;
	Geometry();
	~Geometry();	

	void DestroyScatters();
	void InterpolateVec(Vec vM, Vec vN);
	void CollectVec(Vec vN, Vec vM);
	void Stamp(Vec vN, int ic, int ir, Vec scratchM);
	void TimesI(Vec v, Vec Iv);
	void VecSqMedium(Vec v, Vec vsq, Vec scratchM);

	int Last2(int i);
	void SetJacobian(Mat J, Vec v, int jc, int jr, int jh);

};

class Mode{

	public:

	Vec vpsi;
	std::string name;
	double k[3];
	int ifix, BCPeriod, b[3][2], lasing;
	Mat J;
	KSP ksp; // one ksp per J seems faster

	Mode(Geometry& geo, int ifix_, int b_[3][2], int BCPeriod_, double k_[3]);
	Mode(std::string Name, Geometry& geo, double *Dout); // read constructor
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
void ReadVector(std::ifstream& is, int N, Vec v);


PetscErrorCode MyError(const char* message);
double sqr(const double a);
dcomp sqr(dcomp a);
void View(Mat A, PetscViewer viewer);
void View(Vec x, PetscViewer viewer);

// odd, if I put this in Tools.c, it works for Vec but not Mat. Perhaps something about the Petsc type differences between the two
template<class T> void Output(T A, const char* name, const char* variable_name = ""){

	char filename[PETSC_MAX_PATH_LEN];
	sprintf(filename, "%s%s", name, Output_Suffix.c_str());

	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
	PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject) A, strncmp(variable_name, "", PETSC_MAX_PATH_LEN) ?variable_name : name);
	View(A, viewer);
	PetscViewerDestroy(&viewer);
}


void NewtonSolve(modelist &L, Geometry& geo, Vec v, Vec f, Vec dv);
void ThresholdSearch(double wimag_lo, double wimag_hi, double D_lo, double D_hi, modelist &Lh, Vec vNh, Mode& m, Geometry& geo, Vec f, Vec dv);
double FormJf(modelist& L, Geometry& geo, Vec v, Vec f);


class Vecfun{
// always Nxyzcr()+2!	
	private:
	Vec u;
	int ms, me; // should keep these, since they are used by the val() function
				// every single time.

	protected: 
	double* a;

	public:
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
	double val(int i) const{ return a[i-ns()];}
	void set(int i, double val){ a[i-ns()] = val;}

};

class ComplexVecfun: public Vecfun{

	private:
	Vec v;
	double* b;
	int Nxyzc;
	
	public:
	ComplexVecfun(Vec w, Vec x): Vecfun(w){
		v = x;
		VecGetArray(v, &b);
		
		int N;
		VecGetSize(v, &N);
		Nxyzc = (N-2)/2;
	}
	dcomp val(int i) const{
	
		dcomp z( a[i-ns()], -b[i-ns()] );
		if(i/Nxyzc) z*= ComplexI;
		return z;
	
	}
	void set(int i, dcomp val){
		
		if(i/Nxyzc) val /= ComplexI;
		a[i-ns()] = val.real();
		b[i-ns()] = -val.imag();
	}
	~ComplexVecfun(){ VecRestoreArray(v, &b); }

};


dcomp pmlval(int i, int* N, int* Npml, double* h, int LowerPML, int k);
void AddPlaceholders(Mat J, Geometry &geo);
void AllocateJacobian(Mat J, Geometry& geo);
void AddRowDerivatives(Mat J, Geometry& geo, int ifix, int ih=0);



