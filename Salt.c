#include "headers.h"
#include <algorithm>



PetscErrorCode ReadModes(modelist& L, Geometry& geo){

	char modename[PETSC_MAX_PATH_LEN], 
	     optionin[PETSC_MAX_PATH_LEN] = "-in0",
	     optionout[PETSC_MAX_PATH_LEN] = "-out0";
	for(int i=0; OptionsGet(optionin, modename); i++){

		double D;
		Mode* m = new Mode(modename, geo, &D);

		if(i==0) geo.D = D;
		else if(D != geo.D)
			MyError("The input modes should all be at the same pump strength!");		

		if( !OptionsGet(optionout, modename) )
			MyError("number of -out less than number of -in!");
		m->name = modename;

		L.push_back(m);
		m->Setup(geo);

		sprintf(optionin, "-in%i", i+1);
		sprintf(optionout, "-out%i", i+1);
	}
	return 0;

}



struct yes_lasing{ bool operator() (Mode *m){return m->lasing;}  };
struct not_lasing{ bool operator() (Mode *m){return !m->lasing;}  };
// can put this inside a function too
// btw, this syntax is very interesting. Study this later

int CountLasing(modelist L){ // makes copy
	
	L.remove_if(not_lasing());
	return L.size();
}



void FirstStep(modelist Lh, Mode *m, Geometry& geo, Vec vNh, Vec f, Vec dv, double c){


	PetscPrintf(PETSC_COMM_WORLD, "Taking first step for mode \"%s\"...\n", m->name.c_str() );

  int ih=0;
  FORMODES(Lh, it){ // find ih of m
		if( *it == m) break;
		else ih++;
  }

	if(vNh != m->vpsi){ // update vpsi's from v
		int ih =0;
		FORMODES(Lh, it){
			ScatterRange((*it)->vpsi, vNh, 0, ih*geo.NJ(), geo.NJ() );
			ih++;
		}
	}

  while(1){

	if( LastProcess() ){ // try new c
		VecSetValue(vNh, geo.offset(ih)+geo.Nxyzcr()+1, c, INSERT_VALUES);	
		if( vNh != m->vpsi) VecSetValue(m->vpsi, geo.Nxyzcr()+1, c, INSERT_VALUES);
	}
	Assemble(vNh);
	Assemble(m->vpsi);

	if( FormJf(Lh, geo, vNh, f) < OptionsDouble("-newtonf_tol")) break;
	KSPSolve( (*Lh.begin())->ksp, f, dv);
	PetscPrintf(PETSC_COMM_WORLD, "\n");

	double dc = -GetValue(dv, geo.offset(ih)+geo.Nxyzcr()+1 );

	if( std::abs(dc)/c < 0.5){
		VecAXPY(vNh, -1.0, dv);

		if(vNh != m->vpsi){ // update vpsi's from v
			int ih =0;
			FORMODES(Lh, it){
				ScatterRange(vNh, (*it)->vpsi, ih*geo.NJ(), 0, geo.NJ() );
				ih++;
			}
		}
		// don't NewtonSolve here, that will be done immediately after
		break;
	}else if(c + dc < 0) c *= 0.5;
	else c = 0.5*(c + c+dc);

  }
  
  	PetscPrintf(PETSC_COMM_WORLD, "First step for mode \"%s\" complete!\n", m->name.c_str() );  	


}




void Bundle(modelist &L, Geometry &geo){




	int Nh = L.size(), Nj = 2*geo.Nxyzc()+2;
	if(Nh < 2) MyError("Bundle function is only for multimode!");
	
	Mat J; KSP ksp;
	CreateSquareMatrix( Nh*Nj, 0, &J);
	AllocateJacobian(J, geo);
	
	AddPlaceholders(J, geo);


	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetFromOptions(ksp);
	KSPSetOperators(ksp, J, J, SAME_PRECONDITIONER);
	// TODO: will probably want to merge all of this in with a generalized
	// multimode version of Mode::Setup
	

	for(int i=0; i<SCRATCHNUM; i++){
		Destroy(&geo.vNhscratch[i]);
		MatGetVecs(J, &geo.vNhscratch[i], NULL);
	}


	int ih = 0;

	FORMODES(L, it){
		Mode *m = *it;
		
		Destroy( &m->J); // bundle shares J and v
		m->J = J;
		KSPDestroy(&m->ksp);
		m->ksp = ksp;
		
		geo.MoperatorGeneralBlochFill(J, m->b, m->BCPeriod, m->k, ih);
		AddRowDerivatives(J, geo, m->ifix, ih);
		ih++;
	}	
	
	Assemble(J);
	MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	MatStoreValues(J); 	
	
	
}


bool AtThreshold(Mode *m){ return m->c() == 0.0 && m->lasing;}

int main(int argc, char** argv){ SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); {
	Geometry geo;


	double dD, Dmax; 
	OptionsGet("-dD", &dD);
	OptionsGet("-Dmax", &Dmax);


	
	modelist L;
	ReadModes(L, geo);
	Mode* m = *(L.begin());
	
    Vec f, dv;
    MatGetVecs( m->J, &dv, &f);

	for(; geo.D <= Dmax; geo.D = (geo.D+dD < Dmax? geo.D+dD: Dmax)){

	  modelist Ll = L, Ln = L; // lasing and nonlasing
	  Ll.remove_if(not_lasing() );
	  Ln.remove_if(yes_lasing() );
	  
	  Vec vNh = (*L.begin() )->vpsi, fNh = f, dvNh = dv;

	  if( Ll.size() > 0){ // lasing modes
		  modelist::iterator it = std::find_if(Ll.begin(), Ll.end(), AtThreshold);
		  if( it != Ll.end() && Ll.size() > 1) Bundle(Ll, geo);

		  if(Ll.size() > 1){	 // these vectors will have been properly created in the last block
			vNh = geo.vNhscratch[2];
			fNh = geo.vNhscratch[3];
			dvNh = geo.vNhscratch[4];
		  }

		  if( it != Ll.end() ){
		  		geo.D += 0.5*dD;
		  		FirstStep(Ll, *it, geo, vNh, fNh, dvNh, 1.0);
		  }
		  NewtonSolve(Ll, geo,  vNh, fNh, dvNh);  
	  }

	  FORMODES(Ln, it){
		m = *it;

		modelist Lm;
		Lm.push_back(m); // for non-lasing modes, do one by one

		double wi_old = m->w().imag();
		NewtonSolve(Lm, geo,  m->vpsi, f, dv);
		double wi_new = m->w().imag();

		if(wi_new > -OptionsDouble("-thresholdw_tol") && !m->lasing)
			ThresholdSearch(  wi_old, wi_new, geo.D-dD, geo.D, 
			Ll, vNh, *m, geo, f, dv); // todo: replace with vNh

	  }
	  
	  if(geo.D==Dmax) break;
	}

	FORMODES(L, it){
		(*it)->Write(geo);
		delete *it; // since the pointer was allocated using new
		// also, delete automatically calls destructor
	}

	Destroy(&f);
	Destroy(&dv);
	PetscPrintf(PETSC_COMM_WORLD, "\n");
	PetscPrintf(PETSC_COMM_WORLD, "TODO: a whole bunch of TODOs in Salt.c related to first step of multimode\n");	
	PetscPrintf(PETSC_COMM_WORLD, "future todo: add artificial crashes to enforce all the assumptions I'm making. For example, crash if any file read fails.\n");		
	PetscPrintf(PETSC_COMM_WORLD, "future todo: Make sure all MyError crashes crash all processes the way CHKERRQ does\n");		
	PetscPrintf(PETSC_COMM_WORLD, "TODO: compare multimode calculation speed to single mode with twice the pixels\n");		
}SlepcFinalize();	}
