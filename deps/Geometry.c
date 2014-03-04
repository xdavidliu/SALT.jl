#include "headers.h"

void SetPump(Geometry geo, double D){
	geo->D = D;
}

void InterpolateVec(Geometry geo, Vec vM, Vec vN){
	VecSet(vN, 0.0);

	const double *vals;
	VecGetArrayRead(vM, &vals);
	int ms, me;
	VecGetOwnershipRange(vM, &ms, &me);

	int i;
	for(i=0; i< Nxyzcr(geo); i++){
	
		Point p;
		CreatePoint_i(&p, i, &geo->gN);
		if( projectmedium(&p, &geo->gM, geo->LowerPML) ){
			int ixyz = xyz(&p);
			if( ms <= ixyz && ixyz < me)
				VecSetValue(vN, i, vals[ixyz-ms], ADD_VALUES);
		}
	}

	VecRestoreArrayRead(vM, &vals);
	AssembleVec(vN);
}

void CollectVec(Geometry geo, Vec vN, Vec vM){
	VecSet(vM, 0.0);

	const double *vals;
	VecGetArrayRead(vN, &vals);
	int ns, ne;
	VecGetOwnershipRange(vN, &ns, &ne);
	if( ne > Nxyzcr(geo)-2) ne = Nxyzcr(geo)-2;

	int i;
	for(i=ns; i<ne; i++){
	
		Point p;
		CreatePoint_i(&p, i, &geo->gN);
		if( projectmedium(&p, &geo->gM, geo->LowerPML) )
			VecSetValue(vM, xyz(&p), vals[i-ns], ADD_VALUES);
	}

	VecRestoreArrayRead(vN, &vals);
	AssembleVec(vM);
}

void TimesI(Geometry geo, Vec v, Vec Iv){
	VecSet(Iv, 0.0);
	ScatterRange (v, Iv, Nxyzc(geo), 0, Nxyzc(geo));

	VecScale(Iv, -1.0);
	ScatterRange (v, Iv, 0, Nxyzc(geo), Nxyzc(geo));
}

void VecSqMedium(Geometry geo, Vec v, Vec vsq, Vec scratchM){
	VecCopy(v, vsq);
	VecPointwiseMult(vsq, vsq, vsq);

	CollectVec(geo, vsq, scratchM);
	InterpolateVec(geo, scratchM, vsq);
}

Geometry CreateGeometry(int N[3], double h[3], int Npml[3], int Nc, int LowerPML, double *eps, double *fprof, double wa, double y){
	int i;

        Geometry geo = (Geometry) malloc(sizeof(struct Geometry_s));
	geo->Nc = Nc;
	geo->LowerPML = LowerPML;

	for(i=0; i<3; i++){
		geo->h[i] = h[i];
		geo->Npml[i] = Npml[i];
	}

	CreateGrid(&geo->gN, N, geo->Nc, 2);
	CreateGrid(&geo->gM, N, 1, 1); // 3/3/14: set M = N as per Steven

	CreateVec(2*Nxyzc(geo)+2, &geo->vepspml);

	Vecfun pml;
	CreateVecfun(&pml,geo->vepspml);

	for(i=pml.ns; i<pml.ne; i++){
		Point p;
		CreatePoint_i(&p, i, &geo->gN);
		project(&p, 3);
		dcomp eps_geoal = pmlval(xyzc(&p), N, geo->Npml, geo->h, geo->LowerPML, 0);
		setr(&pml, i, p.ir? cimag(eps_geoal) : creal(eps_geoal) );
	}
	DestroyVecfun(&pml);

	CreateVec(Mxyz(geo), &geo->vMscratch[0]);

	for(i=1; i<SCRATCHNUM; i++) VecDuplicate(geo->vMscratch[0], &geo->vMscratch[i]);

	{ 	double *scratch;
		int ns, ne;
		VecGetOwnershipRange(geo->vMscratch[0], &ns, &ne);
		VecGetArray(geo->vMscratch[0], &scratch);

		for(i=ns; i<ne;i++)
			scratch[i-ns] = eps[i-ns];

		VecRestoreArray(geo->vMscratch[0], &scratch);
	}

	CreateVec(2*Nxyzc(geo)+2, &geo->vH);
	VecDuplicate(geo->vH, &geo->veps);
	VecDuplicate(geo->vH, &geo->vIeps);
	for(i=0; i<SCRATCHNUM; i++) VecDuplicate(geo->vH, &geo->vscratch[i]);
	VecSet(geo->vH, 1.0);

	VecShift(geo->vMscratch[0], -1.0); //hack, for background dielectric
	InterpolateVec(geo, geo->vMscratch[0], geo->vscratch[1]);
	VecShift(geo->vscratch[1], 1.0);

	VecPointwiseMult(geo->veps, geo->vscratch[1], geo->vepspml);
	TimesI(geo, geo->veps, geo->vIeps); // vIeps for convenience only, make sure to update it later if eps ever changes!

	geo->D = 0.0;

	geo->wa = wa;
	geo->y = y;

	VecDuplicate(geo->veps, &geo->vf);

	{ 	double *scratch;
		int ns, ne;
		VecGetOwnershipRange(geo->vMscratch[1], &ns, &ne);
		VecGetArray(geo->vMscratch[1], &scratch);

		for(i=ns; i<ne;i++)
			scratch[i-ns] = fprof[i-ns];

		VecRestoreArray(geo->vMscratch[1], &scratch);
	}

	InterpolateVec(geo, geo->vMscratch[1], geo->vf);

        return geo;
}

void DestroyGeometry(Geometry geo){
    if (geo) {
		int i; 
		for(i=0; i<SCRATCHNUM; i++){
			DestroyVec(&geo->vscratch[i]);
			DestroyVec(&geo->vNhscratch[i]);	
			DestroyVec(&geo->vMscratch[i]);
		}
		DestroyVec(&geo->vH);
		DestroyVec(&geo->veps);
		DestroyVec(&geo->vIeps);
		DestroyVec(&geo->vepspml);
		DestroyVec(&geo->vf);
		free(geo);
    }
}

int Last2(Geometry geo, int i){
	return i%(Nxyzcr(geo)+2) / Nxyzcr(geo);
}

void SetJacobian(Geometry geo, Mat J, Vec v, int jc, int jr, int jh){
// use jc = -1 for last 2 columns, and jc = -2 for blocks (all jc)

	AssembleVec(v);

	int ns, ne, i;
	VecGetOwnershipRange(v, &ns, &ne);
	const double *vals;
	VecGetArrayRead(v, &vals);

	for(i=ns; i<ne; i++){
		if( Last2(geo, i) || vals[i-ns] == 0.0) continue;

		int col, offset = jh*(Nxyzcr(geo)+2),
		ij = i%(Nxyzcr(geo)+2);
	
		Point p;
		CreatePoint_i(&p, ij, &geo->gN);
	
		if(jc == -1) //columns
			col = Nxyzcr(geo)+jr;
		else if(jc == -2) //blocks
			col = jr*Nxyzc(geo) + xyzc(&p);
		else // tensors
			col = jr*Nxyzc(geo) + jc*Nxyz(geo) + xyz(&p);
	
		MatSetValue(J, i, col+offset, vals[i-ns], ADD_VALUES);
	
	}
	VecRestoreArrayRead(v, &vals);
}

// ================= accessor functions for Julia ===========

int GetN(Geometry geo, int i){ return geo->gN.N[i]; }
int GetNpml(Geometry geo, int i){ return geo->Npml[i]; }
double GetCellh(Geometry geo, int i){ return geo->h[i]; }
int GetNc(Geometry geo){ return geo->Nc;}

Vec GetVeps(Geometry geo){ return geo->veps;}
Vec GetVfprof(Geometry geo){ return geo->vf;}
