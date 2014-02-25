#include "headers.h"

void VecSetComplex(Vec vR, Vec vI, int i, int ir, dcomp val, InsertMode addv){

		VecSetValue(vR, i, ir? val.imag() : val.real(), addv );
		VecSetValue(vI, i, ir? val.real() : -val.imag(), addv );
}


void Isolate(Vec v, Grid& gN, int ic, int ir){

	int ns, ne;
	VecGetOwnershipRange(v, &ns, &ne);
	double *a;
	VecGetArray(v, &a);

	for(int i=ns; i<ne && i<xyzcrGrid(&gN); i++){

		Point p(i, gN);
		if(p.ic != ic || p.ir != ir) a[i-ns] = 0.0;
	}
	VecRestoreArray(v, &a);

}


void Geometry::Stamp(Vec vN, int ic, int ir, Vec scratchM){

	Isolate(vN, gN, ic, ir);
	CollectVec(vN, scratchM);
	InterpolateVec(scratchM, vN);

}

void LinearDerivative(Mode& m, Geometry& geo, Vec dfR, Vec dfI, int ih){

	ComplexVecfun eps;
	CreateComplexVecfun(&eps, geo.veps, geo.vIeps);
	Vecfun f(geo.vf), H(geo.vH);

	dcomp mw = m.w(), yw = m.gamma_w(geo);

	for(int i=eps.ns; i<eps.ne; i++){
		dcomp val = sqr(mw) * (valc(&eps, i) + geo.D * yw * f.valr(i) * H.valr(i) );
		VecSetComplex(dfR, dfI, i+offset(&geo, ih), ir(&geo, i), val, INSERT_VALUES);
	}

	DestroyComplexVecfun(&eps);
}

void TensorDerivative(Mode& m, Mode &mj, Geometry& geo, int jc, int jr, Vec df, Vec vpsibra, Vec vIpsi, int ih){


	double mjc = mj.c();
	dcomp mw = m.w(), yw = m.gamma_w(geo), yjw = mj.gamma_w(geo);

	Vecfun f(geo.vf), H(geo.vH);
	ComplexVecfun psi;
	CreateComplexVecfun(&psi, m.vpsi, vIpsi);
	Vecfun psibra(vpsibra);

	for(int i=f.ns(); i<f.ne(); i++){

		if( f.valr(i) == 0.0) continue;			
		dcomp ket_term = -sqr(mw ) * sqr(mjc) * sqr(std::abs(yjw)) * 2.0
			* sqr(H.valr(i) ) * geo.D * f.valr(i) * yw * valc(&psi, i);		
		double val = psibra.valr(i) * (ir(&geo, i)? ket_term.imag() : ket_term.real() );
		
		VecSetValue(df, i+offset(&geo, ih), val, INSERT_VALUES);
	}
	DestroyComplexVecfun(&psi);
}


void ColumnDerivative(Mode* m, Mode* mj, Geometry& geo, Vec dfR, Vec dfI, Vec vIpsi, Vec vpsisq, int ih){
	// vIpsi is for m, vpsisq is for mj	
	// use pointers so can check whether ih = jh

	// purposely don't set df = 0 here to allow multiple ih's
	double mjc = mj->c();
	dcomp mw = m->w(), yw = m->gamma_w(geo),
			mjw = mj->w(), yjw = mj->gamma_w(geo);
	
	ComplexVecfun psi, eps;
	CreateComplexVecfun(&psi,m->vpsi, vIpsi);
	CreateComplexVecfun(&eps,geo.veps, geo.vIeps);

	Vecfun f(geo.vf), H(geo.vH), psisq(vpsisq); 

	for(int i=psi.ns; i<psi.ne; i++){

		dcomp dfdk = 0.0, dfdc = 0.0, 
			DfywHpsi = geo.D * f.valr(i) * yw * H.valr(i) * valc(&psi, i);

		if(m == mj)
			dfdk += ( -sqr(mw)*yw / geo.y +2.0*mw ) * DfywHpsi + 2.0*mw* valc(&eps, i)*valc(&psi, i);

		if(m->lasing && f.valr(i) != 0.0){
			dcomp dHdk_term = -sqr(mjc) * -2.0*(mjw-geo.wa)
			 /sqr(geo.y) * sqr(sqr(std::abs(yjw)));
			dHdk_term *= sqr(mw)*DfywHpsi * H.valr(i) * psisq.valr(i);
			dfdk += dHdk_term;

			dfdc = sqr(mw) * DfywHpsi * H.valr(i);
			dfdc *= (-2.0*mjc)*sqr(std::abs(yjw)) * psisq.valr(i);
		}
		
		if( !m->lasing)
			VecSetComplex(dfR, dfI, i+offset(&geo, ih), ir(&geo, i), dfdk, INSERT_VALUES);
		else{
			VecSetValue(dfR, i+offset(&geo, ih), ir(&geo, i)? dfdk.imag() : dfdk.real(), INSERT_VALUES );
			VecSetValue(dfI, i+offset(&geo, ih), ir(&geo, i)? dfdc.imag() : dfdc.real(), INSERT_VALUES );
		}
	}

	DestroyComplexVecfun(&eps);
	DestroyComplexVecfun(&psi);

}

void ComputeGain(Geometry& geo, modelist& L){

  VecSet(geo.vH, 0.0);

  Vecfun H(geo.vH);
  FORMODES(L, it){
	Mode *m = *it;
	dcomp yw = m->gamma_w(geo);
	double mc = m->c();

	// do not change this from vscratch[3], or the hack below for single mode Column derivative will fail!
	geo.VecSqMedium(m->vpsi, geo.vscratch[3], geo.vMscratch[0]);

	Vecfun psisq(geo.vscratch[3]);
	for(int i=H.ns(); i<H.ne(); i++)
		H.setr(i, H.valr(i) + sqr(std::abs(yw)) *sqr(mc) * psisq.valr(i) ) ;
   }

   for(int i=H.ns(); i<H.ne(); i++)
		H.setr(i, 1.0 / (1.0 + H.valr(i) ) );
}

double FormJf(modelist& L, Geometry& geo, Vec v, Vec f){

	Mode *m = *(L.begin());
	bool lasing = m->lasing;
	Mat J = m->J; // for multimode, all m share same J

	if(lasing)
		ComputeGain(geo, L); // do this before naming scratch vectors!

	// ================== name scratch vectors ================== //

	Vec vpsisq = geo.vscratch[3], // only form this later if needed
		vIpsi = geo.vscratch[2];

	Vec dfR, dfI;
	if(L.size() == 1){
		dfR = geo.vscratch[0];
		dfI = geo.vscratch[1];
	}else{
		dfR = geo.vNhscratch[0];
		dfI = geo.vNhscratch[1];
	}


	// =========== linear J to compute residual ========= //
	MatRetrieveValues(J);

	int ih, jh;

	ih = 0;
	FORMODES(L, it){
		m = *it;
		VecSet(dfR, 0.0);		
		VecSet(dfI, 0.0);
	
		LinearDerivative(*m, geo, dfR, dfI, ih);

	  	geo.SetJacobian(J, dfR, -2, 0, ih);
		geo.SetJacobian(J, dfI, -2, 1, ih); 

		ih++;
	}

	// row derivatives already added in add placeholders!

	AssembleMat(J);

	MatMult(J, v, f);

	for(int kh = 0; kh<L.size(); kh++) for(int ir=0; ir<2; ir++)
		VecSetValue(f, kh*NJ(&geo) + Nxyzcr(&geo)+ir, 0.0, INSERT_VALUES);

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
			geo.VecSqMedium(mj->vpsi, vpsisq, geo.vMscratch[0]);

		ih = 0;
		FORMODES(L, it){
			Mode *mi = *it;

			geo.TimesI(mi->vpsi, vIpsi);
			ColumnDerivative(mi, mj, geo, dfR, dfI, vIpsi, vpsisq, ih);

			ih++;
		}

		geo.SetJacobian(J, dfR, -1, 0, jh);
		geo.SetJacobian(J, dfI, -1, 1, jh);


		jh++;
	}

	//================ tensor derivatives ================
	
  if(lasing){
	Vec vpsibra = vpsisq; vpsisq = 0;

	jh = 0;
	FORMODES(L, jt){
		Mode *mj = *jt;

		for(int jr=0; jr<2; jr++) for(int jc=0; jc< geo.gN.Nc; jc++){

			VecCopy(mj->vpsi, vpsibra);
			geo.Stamp(vpsibra, jc, jr, geo.vMscratch[0]);

			VecSet(dfR, 0.0);
			ih = 0;
			FORMODES(L, it){
				Mode *mi = *it;
				geo.TimesI(mi->vpsi, vIpsi);
	
				TensorDerivative(*mi, *mj, geo, jc, jr, dfR, vpsibra, vIpsi, ih);
				ih++;
			}

			geo.SetJacobian(J, dfR, jc, jr, jh);
		}

		jh++;
	}
  }
	AssembleMat(J);


	return fnorm;	

}
