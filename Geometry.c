#include "headers.h"





void Geometry::InterpolateVec(Vec vM, Vec vN){

	VecSet(vN, 0.0);

	const double *vals;
	VecGetArrayRead(vM, &vals);	
	int ms, me;
	VecGetOwnershipRange(vM, &ms, &me);

	for(int i=0; i< Nxyzcr(this); i++){
		
		Point p(i, gN);
		if( projectmedium(&p, gM, LowerPML) )
			VecSetValue(vN, i, vals[xyz(&p)-ms], ADD_VALUES);
	}

	VecRestoreArrayRead(vM, &vals);	
	AssembleVec(vN);

}

void Geometry::CollectVec(Vec vN, Vec vM){

	VecSet(vM, 0.0);

	const double *vals;
	VecGetArrayRead(vN, &vals);
	int ns, ne;
	VecGetOwnershipRange(vN, &ns, &ne);
	if( ne > Nxyzcr(this)-2) ne = Nxyzcr(this)-2;
	
	for(int i=ns; i<ne; i++){
		
		Point p(i, gN);
		if( projectmedium(&p, gM, LowerPML) )
			VecSetValue(vM, xyz(&p), vals[i-ns], ADD_VALUES);
	}
	
	VecRestoreArrayRead(vN, &vals);
	AssembleVec(vM);
}


void Geometry::TimesI(Vec v, Vec Iv){

	VecSet(Iv, 0.0);
	ScatterRange (v, Iv, Nxyzc(this), 0, Nxyzc(this));
	
	
	VecScale(Iv, -1.0);
	ScatterRange (v, Iv, 0, Nxyzc(this), Nxyzc(this));

}

void Geometry::VecSqMedium(Vec v, Vec vsq, Vec scratchM){
		VecCopy(v, vsq);
		VecPointwiseMult(vsq, vsq, vsq);

		CollectVec(vsq, scratchM);
		InterpolateVec(scratchM, vsq);
}


Geometry::Geometry(){


	int N[3], M[3];

	OptionsXYZInt("-N", N);
	OptionsXYZInt("-M", M);
	OptionsXYZInt("-Npml", Npml);
	OptionsXYZDouble("-h", h);

	OptionsGetInt("-Nc", &Nc);
	OptionsGetInt("-LowerPML", &LowerPML);
	gN = Grid(N, Nc, 2);
	gM = Grid(M, 1, 1);


	CreateVec(2*Nxyzc(this)+2, &vepspml);



	Vecfun pml(vepspml);

	for( int i=pml.ns; i<pml.ne; i++){
		Point p(i, gN);
		project(&p, 3);
		dcomp eps_geoal = pmlval(xyzc(&p), N, Npml, h, LowerPML, 0);	
		setr(&pml, i, p.ir? eps_geoal.imag() : eps_geoal.real() );
	}
	DestroyVecfun(&pml);
	
	CreateVec(Mxyz(this), &vMscratch[0]);

	for(int i=1; i<SCRATCHNUM; i++) VecDuplicate(vMscratch[0], &vMscratch[i]);
	char file[PETSC_MAX_PATH_LEN];
	OptionsGetString("-epsfile", file);

	FILE *fp;
	
	fp = fopen(file, "r");
	if(fp==NULL){
		char message[PETSC_MAX_PATH_LEN];
		sprintf(message, "failed to read %s", file);
		MyError(message);
	}
	ReadVectorC(fp, Mxyz(this), vMscratch[0]);
	fclose(fp);
	
	CreateVec(2*Nxyzc(this)+2, &vH);
	VecDuplicate(vH, &veps);
	VecDuplicate(vH, &vIeps);
	for(int i=0; i<SCRATCHNUM; i++) VecDuplicate(vH, &vscratch[i]);
	VecSet(vH, 1.0);	


	VecShift(vMscratch[0], -1.0); //hack, for background dielectric
	InterpolateVec(vMscratch[0], vscratch[1]);
	VecShift(vscratch[1], 1.0);


	VecPointwiseMult(veps, vscratch[1], vepspml);
	TimesI(veps, vIeps); // vIeps for convenience only, make sure to update it later if eps ever changes!



	D = 0.0;
	OptionsGetDouble("-wa", &wa);
	OptionsGetDouble("-gamma", &y);


	VecDuplicate(veps, &vf);

	OptionsGetString("-fproffile", file);

	fp = fopen(file, "r");	
	if(fp==NULL){
		char message[PETSC_MAX_PATH_LEN];
		sprintf(message, "failed to read %s", file);
		MyError(message);
	}	
	ReadVectorC(fp, Mxyz(this), vMscratch[1]);
	fclose(fp);	
	
	InterpolateVec(vMscratch[1], vf);



 
}



Geometry::~Geometry(){

	for(int i=0; i<SCRATCHNUM; i++){
		DestroyVec(&vscratch[i]);
		DestroyVec(&vNhscratch[i]);		
		DestroyVec(&vMscratch[i]);
	}
	DestroyVec(&vH);
	DestroyVec(&veps);
	DestroyVec(&vIeps);
	DestroyVec(&vepspml);
	DestroyVec(&vf);

}

int Geometry::Last2(int i){
	return i%(Nxyzcr(this)+2) / Nxyzcr(this);
}


void Geometry::SetJacobian(Mat J, Vec v, int jc, int jr, int jh){
// use jc = -1 for last 2 columns, and jc = -2 for blocks (all jc)

	AssembleVec(v);

	int ns, ne;
	VecGetOwnershipRange(v, &ns, &ne);
	const double *vals;
	VecGetArrayRead(v, &vals);
	
	for(int i=ns; i<ne; i++){
		if( Last2(i) || vals[i-ns] == 0.0) continue;
	
		int col, offset = jh*(Nxyzcr(this)+2),
		ij = i%(Nxyzcr(this)+2);
		
		Point p(ij, gN);
		
		if(jc == -1) //columns
			col = Nxyzcr(this)+jr;
		else if(jc == -2) //blocks
			col = jr*Nxyzc(this) + xyzc(&p);
		else // tensors
			col = jr*Nxyzc(this) + jc*Nxyz(this) + xyz(&p);
		
		MatSetValue(J, i, col+offset, vals[i-ns], ADD_VALUES);
		
	}
	VecRestoreArrayRead(v, &vals);	

}



