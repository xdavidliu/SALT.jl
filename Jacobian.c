#include "headers.h"

void VecSetComplex(Vec vR, Vec vI, int i, int ir, dcomp val, InsertMode addv){

		VecSetValue(vR, i, ir? val.imag() : val.real(), addv );
		VecSetValue(vI, i, ir? val.real() : -val.imag(), addv );
}


void Isolate(Vec v, Grid& gN, int ic, int ir){

	int ns, ne, i;
	VecGetOwnershipRange(v, &ns, &ne);
	double *a;
	VecGetArray(v, &a);

	for(i=ns; i<ne && i<xyzcrGrid(&gN); i++){

		Point p;
		CreatePoint_i(&p, i, &gN);
		if(p.ic != ic || p.ir != ir) a[i-ns] = 0.0;
	}
	VecRestoreArray(v, &a);

}


void Stamp(Geometry *geo, Vec vN, int ic, int ir, Vec scratchM){

	Isolate(vN, geo->gN, ic, ir);
	CollectVec(geo, vN, scratchM);
	InterpolateVec(geo, scratchM, vN);

}

void LinearDerivative(Mode *m, Geometry *geo, Vec dfR, Vec dfI, int ih){

	Complexfun eps;
	CreateComplexfun(&eps, geo->veps, geo->vIeps);
	Vecfun f, H;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);

	dcomp mw = getw(m), yw = gamma_w(m, geo);

	int i;
	for(i=eps.ns; i<eps.ne; i++){
		dcomp val = sqr(mw) * (valc(&eps, i) + geo->D * yw * valr(&f, i) * valr(&H, i) );
		VecSetComplex(dfR, dfI, i+offset(geo, ih), ir(geo, i), val, INSERT_VALUES);
	}
	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyComplexfun(&eps);
}

void TensorDerivative(Mode *m, Mode *mj, Geometry *geo, int jc, int jr, Vec df, Vec vpsibra, Vec vIpsi, int ih){


	double mjc = getc(mj);
	dcomp mw = getw(m), yw = gamma_w(m, geo), yjw = gamma_w(mj, geo);

	Vecfun f, H, psibra;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&psibra, vpsibra);

	Complexfun psi;
	CreateComplexfun(&psi, m->vpsi, vIpsi);

	int i;
	for(i=f.ns; i<f.ne; i++){

		if( valr(&f, i) == 0.0) continue;			
		dcomp ket_term = -sqr(mw ) * sqr(mjc) * sqr(std::abs(yjw)) * 2.0
			* sqr(valr(&H, i) ) * geo->D * valr(&f, i) * yw * valc(&psi, i);		
		double val = valr(&psibra, i) * (ir(geo, i)? ket_term.imag() : ket_term.real() );
		
		VecSetValue(df, i+offset(geo, ih), val, INSERT_VALUES);
	}
	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyVecfun(&psibra);
	DestroyComplexfun(&psi);
}


void ColumnDerivative(Mode* m, Mode* mj, Geometry *geo, Vec dfR, Vec dfI, Vec vIpsi, Vec vpsisq, int ih){
	// vIpsi is for m, vpsisq is for mj	
	// use pointers so can check whether ih = jh

	// purposely don't set df = 0 here to allow multiple ih's
	double mjc = getc(mj);
	dcomp mw = getw(m), yw = gamma_w(m, geo),
			mjw = getw(mj), yjw = gamma_w(mj, geo);
	
	Complexfun psi, eps;
	CreateComplexfun(&psi,m->vpsi, vIpsi);
	CreateComplexfun(&eps,geo->veps, geo->vIeps);

	Vecfun f,H, psisq;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&psisq, vpsisq);

	int i;
	for(i=psi.ns; i<psi.ne; i++){

		dcomp dfdk = 0.0, dfdc = 0.0, 
			DfywHpsi = geo->D * valr(&f, i) * yw * valr(&H, i) * valc(&psi, i);

		if(m == mj)
			dfdk += ( -sqr(mw)*yw / geo->y +2.0*mw ) * DfywHpsi + 2.0*mw* valc(&eps, i)*valc(&psi, i);

		if(m->lasing && valr(&f, i) != 0.0){
			dcomp dHdk_term = -sqr(mjc) * -2.0*(mjw-geo->wa)
			 /sqr(geo->y) * sqr(sqr(std::abs(yjw)));
			dHdk_term *= sqr(mw)*DfywHpsi * valr(&H, i) * valr(&psisq, i);
			dfdk += dHdk_term;

			dfdc = sqr(mw) * DfywHpsi * valr(&H, i);
			dfdc *= (-2.0*mjc)*sqr(std::abs(yjw)) * valr(&psisq, i);
		}
		
		if( !m->lasing)
			VecSetComplex(dfR, dfI, i+offset(geo, ih), ir(geo, i), dfdk, INSERT_VALUES);
		else{
			VecSetValue(dfR, i+offset(geo, ih), ir(geo, i)? dfdk.imag() : dfdk.real(), INSERT_VALUES );
			VecSetValue(dfI, i+offset(geo, ih), ir(geo, i)? dfdc.imag() : dfdc.real(), INSERT_VALUES );
		}
	}

	DestroyComplexfun(&eps);
	DestroyComplexfun(&psi);
	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyVecfun(&psisq);
}

void ComputeGain(Geometry *geo, modelist& L){

	VecSet(geo->vH, 0.0);

	Vecfun H;
	CreateVecfun(&H, geo->vH);
	int i;
	FORMODES(L, it){
		Mode *m = *it;
		dcomp yw = gamma_w(m, geo);
		double mc = getc(m);

		// do not change this from vscratch[3], or the hack below for single mode Column derivative will fail!
		VecSqMedium(geo, m->vpsi, geo->vscratch[3], geo->vMscratch[0]);

		Vecfun psisq;
		CreateVecfun(&psisq ,geo->vscratch[3]);

	
		for(i=H.ns; i<H.ne; i++)
			setr(&H, i, valr(&H, i) + sqr(std::abs(yw)) *sqr(mc) * valr(&psisq, i) ) ;
	}

	for(i=H.ns; i<H.ne; i++)
		setr(&H, i, 1.0 / (1.0 + valr(&H, i) ) );

	DestroyVecfun(&H);
}

double FormJf(modelist& L, Geometry *geo, Vec v, Vec f){

	Mode *m = *(L.begin());
	bool lasing = m->lasing;
	Mat J = m->J; // for multimode, all m share same J

	if(lasing)
		ComputeGain(geo, L); // do this before naming scratch vectors!

	// ================== name scratch vectors ================== //

	Vec vpsisq = geo->vscratch[3], // only form this later if needed
		vIpsi = geo->vscratch[2];

	Vec dfR, dfI;
	if(L.size() == 1){
		dfR = geo->vscratch[0];
		dfI = geo->vscratch[1];
	}else{
		dfR = geo->vNhscratch[0];
		dfI = geo->vNhscratch[1];
	}


	// =========== linear J to compute residual ========= //
	MatRetrieveValues(J);

	int ih, jh, kh, ir, jr, jc;

	ih = 0;
	FORMODES(L, it){
		m = *it;
		VecSet(dfR, 0.0);		
		VecSet(dfI, 0.0);
	
		LinearDerivative(m, geo, dfR, dfI, ih);

	  	SetJacobian(geo, J, dfR, -2, 0, ih);
		SetJacobian(geo, J, dfI, -2, 1, ih); 

		ih++;
	}

	// row derivatives already added in add placeholders!

	AssembleMat(J);

	MatMult(J, v, f);

	for(kh = 0; kh<L.size(); kh++) for(ir=0; ir<2; ir++)
		VecSetValue(f, kh*NJ(geo) + Nxyzcr(geo)+ir, 0.0, INSERT_VALUES);

	AssembleVec(f);

	double fnorm;
	VecNorm(f, NORM_2, &fnorm);
	
	static int printnewton = OptionsInt("-printnewton");
	if( printnewton ) PetscPrintf(PETSC_COMM_WORLD, "|f| = %.0e;", fnorm);
	// no \n here to make room for timing printf statement immediately afterwards


	if(fnorm < OptionsDouble("-newtonf_tol") )
		return fnorm;   		// TODO: deleted old integral routine. Write new one here.

	// =============== column derivatives ====================


	jh = 0;
	FORMODES(L, jt){
		Mode *mj = *jt;

		VecSet(dfR, 0.0);
		VecSet(dfI, 0.0);

		if(L.size() > 1)  // hack: only recompute vpsisq if ComputeGain didn't already do it, i.e. for multimode
			VecSqMedium(geo, mj->vpsi, vpsisq, geo->vMscratch[0]);

		ih = 0;
		FORMODES(L, it){
			Mode *mi = *it;

			TimesI(geo, mi->vpsi, vIpsi);
			ColumnDerivative(mi, mj, geo, dfR, dfI, vIpsi, vpsisq, ih);

			ih++;
		}

		SetJacobian(geo, J, dfR, -1, 0, jh);
		SetJacobian(geo, J, dfI, -1, 1, jh);


		jh++;
	}

	//================ tensor derivatives ================
	
  if(lasing){
	Vec vpsibra = vpsisq; vpsisq = 0;

	jh = 0;
	FORMODES(L, jt){
		Mode *mj = *jt;

		for(jr=0; jr<2; jr++) for(jc=0; jc< geo->gN.Nc; jc++){

			VecCopy(mj->vpsi, vpsibra);
			Stamp(geo, vpsibra, jc, jr, geo->vMscratch[0]);

			VecSet(dfR, 0.0);
			ih = 0;
			FORMODES(L, it){
				Mode *mi = *it;
				TimesI(geo, mi->vpsi, vIpsi);
	
				TensorDerivative(mi, mj, geo, jc, jr, dfR, vpsibra, vIpsi, ih);
				ih++;
			}

			SetJacobian(geo, J, dfR, jc, jr, jh);
		}

		jh++;
	}
  }
	AssembleMat(J);


	return fnorm;	

}
