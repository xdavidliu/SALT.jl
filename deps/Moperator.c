#include "headers.h"



int periodic(int ic, int BCPeriod){
	return ic == (BCPeriod-1) || BCPeriod == 4 || (BCPeriod<0 && ic !=-(BCPeriod+1));
}


dcomp zval(int i, const double* f, int Nxyzc){
	i %= Nxyzc;
	return f[i] + ComplexI*f[i + Nxyzc];
}

int cyclic(Point *P, int ic, int* N){

	return (P->ic==ic)*N[1]*N[2] + (P->ic==(ic+1)%3 )*N[2] + (P->ic==(ic+2)%3 );
} 

void MoperatorGeneralBlochFill(Geometry geo, Mat A, int b[3][2], int DimPeriod, double k[3], int ih){


	int N[3], i, j, l, ib, w, jr, ic;
	for(i=0; i<3; i++) N[i] = geo->gN.N[i];

	double blochbc[3];
	for(i=0; i<3; i++) blochbc[i] = k[i]*N[i]*geo->h[i];

	int NC = 3, offset = ih*(Nxyzcr(geo)+2);
	int ns, ne;
	double hh;

	int bc[3][2][3]; /* bc[x/y/z direction][lo/hi][Ex/Ey/Ez] */

	dcomp val, magicnum, mucp[2], mulcp[2];
	dcomp cidu_phase, cpidu_phase[2], cpidl_phase[2];

	Grid gc, gC;
	CreateGrid(&gc, N, geo->Nc, 2);
	CreateGrid(&gC, N, 3, 2); 

 
  /* set up b ... */
 
	for(ic=0; ic<3; ic++) for(j=0; j<2; j++) for(l=0; l<3; l++)
		bc[ic][j][l] =  b[ic][j]*( l==ic ? -1 :1);

	
	MatGetOwnershipRange(A, &ns, &ne);

	int itrue;
for (itrue = ns; itrue < ne && itrue < 2*Nxyzc(geo); ++itrue) {

	Point p;
	CreatePoint_i(&p, itrue, &gc); 

	project(&p, 3); int i = xyzcr(&p);
	int cp[2], icp[2], cidu, cpidu[2],cpidl[2], cid, cpid[2];
	for(j=0; j<2;j++){

		cp[j] = (p.ic+1+j) % NC;
		icp[j] = i + (cp[j]-p.ic ) *Nxyz(geo);
		cpidu[j] = cyclic(&p, 2-j, N);
		cpidl[j] = cyclic(&p, 2-j, N);
		cpid[j] = cyclic(&p, 2-j, N);

		cpidu_phase[j] = 1.0;
		cpidl_phase[j] = 1.0;
	}	  
	cidu = cyclic(&p, 0, N);
	cid = cyclic(&p, 0, N);

	cidu_phase = 1.0;
    
   for(jr=0; jr<2; jr++) { /* column real/imag parts */
   	int jrd =  (jr-p.ir)*NC*Nxyz(geo);            
   	magicnum = (p.ir==jr)*1.0 + (p.ir<jr)*1.0*ComplexI - (p.ir>jr)*1.0* ComplexI; 

//=====================================================================
	Point prow;
	CreatePoint_i(&prow, i, &gC); 
	project(&prow, geo->Nc);
for(ib=0; ib<2; ib++){


	if(p.ix[p.ic] == N[p.ic]-1){
		int per = periodic(p.ic, DimPeriod );
		cidu = per ? (1-N[p.ic])*cid : 0;
		cidu_phase = per? cexp(ComplexI*blochbc[p.ic]) : bc[p.ic][1][cp[ib]];
	}

	if(p.ix[cp[ib]] == 0){
		int per = periodic(cp[ib], DimPeriod );
		cpidl[ib] = per ? (1-N[cp[ib]])*cpid[ib] : 0;
		cpidl_phase[ib] = per ? cexp(-ComplexI*blochbc[cp[ib]]) : bc[cp[ib]][0][cp[ib]];
	}

        mucp[1-ib] = pmlval(icp[1-ib], N, geo->Npml, geo->h, geo->LowerPML, 1);
      	mulcp[1-ib] = pmlval(icp[1-ib]-cpidl[ib], N, geo->Npml, geo->h, geo->LowerPML, 1);


	double c[4];
        hh = geo->h[p.ic]*geo->h[cp[ib]];
	val = mucp[1-ib] * magicnum /hh; c[1] = creal(val);
	val *= cidu_phase; c[0] = -creal(val);
	val = cpidl_phase[ib] * mulcp[1-ib] * magicnum/hh; c[3] = -creal(val);
	val *= -cidu_phase; c[2] = -creal(val);
      

	int dcol[4];
	dcol[0] = cidu; 
	dcol[1] = 0;
	dcol[2] = cidu-cpidl[ib];
	dcol[3] = -cpidl[ib];



	for(w=0;w<4;w++){
	Point pcol;
	CreatePoint_i(&pcol, icp[ib] + jrd+dcol[w], &gC );

	project(&pcol, geo->Nc);	
	if(pcol.ic!=-1) MatSetValue(A, xyzcr(&prow)+offset, xyzcr(&pcol)+offset, c[w], ADD_VALUES);
	}


	if(p.ix[cp[ib]] == N[cp[ib]]-1){
		int per = periodic(cp[ib], DimPeriod );
		cpidu[ib] = per ? (1-N[cp[ib]])*cpid[ib] : 0;
		cpidu_phase[ib] = per? cexp(ComplexI*blochbc[cp[ib]]) : bc[cp[ib]][1][p.ic];
	}

	if(p.ix[cp[ib]] == 0){
		int per = periodic(cp[ib], DimPeriod );
		cpidl[ib] = per ? (1-N[cp[ib]])*cpid[ib] : -cpidu[ib];
		cpidl_phase[ib] = per? cexp(-ComplexI*blochbc[cp[ib]]) : bc[cp[ib]][0][p.ic];
	}   

     
	hh =  geo->h[cp[ib]]*geo->h[cp[ib]];
	val = -(cpidu_phase[ib] * mucp[1-ib] * magicnum)/hh; c[0] = -creal(val);
	val = +( (mucp[1-ib] + mulcp[1-ib]) * magicnum)/hh;c[1] = -creal(val);
	val = -(cpidl_phase[ib] * mulcp[1-ib] * magicnum)/hh;c[2] = -creal(val);
     
 
	dcol[0] = cpidu[ib]; 
	dcol[1] = 0;
	dcol[2] = -cpidl[ib];

	for(w=0;w<3;w++){
	Point pcol;
	CreatePoint_i(&pcol, i + jrd+dcol[w], &gC );
	project(&pcol, geo->Nc);
	if(pcol.ic!=-1) MatSetValue(A, xyzcr(&prow)+offset, xyzcr(&pcol)+offset, c[w], ADD_VALUES);
	}

}


    }
  }
}


