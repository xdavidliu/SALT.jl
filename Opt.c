#include "headers.h"
#include <nlopt.hpp>

typedef std::vector<double> vecd;


class Plan{

	Geometry geo;
	Vec dv, f;

	
	double dD;
	
	public:
	Mode* m;
	Plan(){
	
		OptionsGet("-dD", &dD);
	
		char name[PETSC_MAX_PATH_LEN];
		OptionsGet("-optin", name);
		m = new Mode(name, geo, &geo.D);
		OptionsGet("-optout", name);
		m->name = name;
		
		double val;
		GetLast2(m->vpsi, &val, NULL);
		SetLast2(m->vpsi, val, -1.0e-14);	// pretend it's still below threshold.	
		m->lasing = 0;
		
		m->Setup(geo); 
		MatGetVecs(m->J, &dv, &f);
		

	}


	void updatewa(double wa_new){ geo.wa = wa_new; }
	double getwa(){ return geo.wa;}
	double threshold(){ return geo.D; }
	void write(){ m->Write(geo); }

	void Solve(){
	
		modelist L;
		L.push_back(m);	

		NewtonSolve(L, geo,  m->vpsi, f, dv);
		geo.D += dD;  // Threshold search allows both points on same side
	
				
		double wi_old = m->w().imag();
		NewtonSolve(L, geo, m->vpsi, f, dv);
		double wi_new = m->w().imag();		
	
		modelist empty;
		ThresholdSearch(
			wi_old, wi_new, geo.D-dD, geo.D, empty, m->vpsi, *m, 
			geo, f, dv);	

		double val;
		GetLast2(m->vpsi, &val, NULL);
		SetLast2(m->vpsi, val, -1.0e-14);	// pretend it's still below threshold.				
		m->lasing = 0;
	}

	~Plan(){
		
		Destroy(&f);
		Destroy(&dv);
		
		delete m;
	}

};



double f(const vecd& x, vecd& grad, void* data){

	PetscPrintf(PETSC_COMM_WORLD, "Evaluating at wa = %g\n", x[0]);
	Plan* p = (Plan *) data;
	p->updatewa(x[0]);
	
	p->Solve();
	
	return p->threshold();
}


int main(int argc, char** argv){ SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); {

	
	Plan p;
	vecd x(1); x[0] = p.getwa();;	

	nlopt::opt opt(nlopt::LN_COBYLA, 1);
	opt.set_min_objective(f, &p);
	opt.set_xtol_rel(1e-2);
	opt.set_initial_step(1.0);
	
	double minf;
	opt.optimize(x, minf);	
	p.write();

}	SlepcFinalize();}