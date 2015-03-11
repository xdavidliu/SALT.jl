#include "salt.h"

void VecSetComplex(Vec vR, Vec vI, int i, int ir, dcomp val, InsertMode addv){
	VecSetValue(vR, i, ir? cimag(val) : creal(val), addv );
	VecSetValue(vI, i, ir? creal(val) : -cimag(val), addv );
}

void Isolate(Vec v, Grid *gN, int ic, int ir){
	int ns, ne, i;
	VecGetOwnershipRange(v, &ns, &ne);
	double *a;
	VecGetArray(v, &a);

	for(i=ns; i<ne && i<xyzcrGrid(gN); i++){
		Point p;
		CreatePoint_i(&p, i, gN);
		if(p.ic != ic || p.ir != ir) a[i-ns] = 0.0;
	}
	VecRestoreArray(v, &a);
}

void Stamp(Geometry geo, Vec vN, int ic, int ir, Vec scratchM){
	Isolate(vN, &geo->gN, ic, ir);
	CollectVec(geo, vN, scratchM);
	InterpolateVec(geo, scratchM, vN);
}

void LinearDerivative(Mode m, Geometry geo, Vec dfR, Vec dfI, int ih){

	Complexfun eps;
	CreateComplexfun(&eps, geo->veps, geo->vIeps);
	Vecfun f, H;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);

	dcomp mw = get_w(m), yw = gamma_w(m, geo);

	int i;
	for(i=eps.ns; i<eps.ne; i++){
		dcomp val = csqr(mw) * (valc(&eps, i) + geo->D * yw * valr(&f, i) * valr(&H, i) );
		VecSetComplex(dfR, dfI, i+offset(geo, ih), ir(geo, i), val, INSERT_VALUES);
		// df is assembled in SetJacobian
	}

	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyComplexfun(&eps);

}

void TensorDerivative(Mode m, Mode mj, Geometry geo, Vec df, Vec vpsibra, Vec vIpsi, int ih){
	double mjc = get_c(mj);
	dcomp mw = get_w(m), yw = gamma_w(m, geo) ;

	Vecfun f, H, psibra;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&psibra, vpsibra);

	Complexfun psi;
	CreateComplexfun(&psi, m->vpsi, vIpsi);

	int i;
	for(i=f.ns; i<f.ne; i++){
		if( valr(&f, i) == 0.0) continue;		
		dcomp ket_term = -csqr(mw ) * sqr(mjc) * 2.0
			* sqr(valr(&H, i) ) * geo->D * valr(&f, i) * yw * valc(&psi, i);	
		double val = valr(&psibra, i) * (ir(geo, i)? cimag(ket_term) : creal(ket_term) );
	
		VecSetValue(df, i+offset(geo, ih), val, INSERT_VALUES);
		// df is assembled in SetJacobian
	}
	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyVecfun(&psibra);
	DestroyComplexfun(&psi);
}

void ColumnDerivative(Mode m, Mode mj, Geometry geo, Vec dfR, Vec dfI, Vec vIpsi, Vec vpsisq, int ih){
	// vIpsi is for m, vpsisq is for mj
	// use pointers so can check whether ih = jh


	// purposely don't set df = 0 here to allow multiple ih's
	double mjc = get_c(mj);
	dcomp mw = get_w(m), yw = gamma_w(m, geo);

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
			dfdk += ( -csqr(mw)*yw / geo->y +2.0*mw ) * DfywHpsi + 2.0*mw* valc(&eps, i)*valc(&psi, i);
		// note: adding dcomp to a double ignores the imaginary part

		if(m->lasing && valr(&f, i) != 0.0){
		
			// dHdk removed; field simply rescaled -DL 6/15/14

			dfdc = csqr(mw) * DfywHpsi * valr(&H, i);
			dfdc *= (-2.0*mjc)*valr(&psisq, i);
		}
	
		if( !m->lasing)
			VecSetComplex(dfR, dfI, i+offset(geo, ih), ir(geo, i), dfdk, INSERT_VALUES);
		else{
			VecSetValue(dfR, i+offset(geo, ih), ir(geo, i)? cimag(dfdk) : creal(dfdk), INSERT_VALUES );
			VecSetValue(dfI, i+offset(geo, ih), ir(geo, i)? cimag(dfdc) : creal(dfdc), INSERT_VALUES );
		// df is assembled in SetJacobian
		}
	}


	DestroyComplexfun(&eps);
	DestroyComplexfun(&psi);
	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyVecfun(&psisq);
}

void ComputeGain(Geometry geo, Mode *ms, int Nh){	

	VecSet(geo->vH, 0.0);
	Vecfun H, f;
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&f, geo->vf);
	int i, ih;
	for(ih=0; ih<Nh; ih++){
		Mode m = ms[ih];
		double mc = get_c(m);

		// do not change this from vscratch[3], or the hack below for single mode Column derivative will fail!
		VecDotMedium(geo, m->vpsi, m->vpsi, geo->vscratch[3], geo->vMscratch[0]);

		Vecfun psisq;
		CreateVecfun(&psisq ,geo->vscratch[3]);
		for(i=H.ns; i<H.ne; i++){
			if(valr(&f, i) == 0.0) continue;
			setr(&H, i, valr(&H, i) + sqr(mc) * valr(&psisq, i) ) ;
		}
		DestroyVecfun(&psisq);
	}
	
	if(geo->interference != 0.0 && Nh == 2){
		// does not affect single mode case
		VecDotMedium(geo, ms[0]->vpsi, ms[1]->vpsi, geo->vscratch[3], geo->vMscratch[0]);

		Vec Ipsi = geo->vscratch[5];
		TimesI( geo, ms[0]->vpsi, Ipsi);
		VecDotMedium(geo, ms[1]->vpsi, Ipsi, geo->vscratch[6], geo->vMscratch[0]);

		// 2 c1 c2 Re[ exp(i thet) psi1* x psi2) ]
		// term in square bracket is (cos thet + i sin thet ) x 
		// ( E1R . E2R + E1I . E2I ) + i ( E1R . E2I - E1I . E2R )
		
		// vscratch[3] and vscratch[6] are the real and imaginary parts of this last line
		double costh = cos(geo->interference), sinth = sin(geo->interference);
		VecScale(geo->vscratch[6], -sinth);
		VecAXPY( geo->vscratch[6], costh, geo->vscratch[3]);
		// now vscratch[6] = Re[ ... ]

		double mc[2] = {get_c(ms[0]), get_c(ms[1]) };
		Vecfun psi_int;
		CreateVecfun(&psi_int ,geo->vscratch[6]);

		for(i=H.ns; i<H.ne; i++){
			if(valr(&f, i) == 0.0) continue;
			setr(&H, i, valr(&H, i) + 2.0*mc[0]*mc[1] * valr(&psi_int, i) ) ;
		}
		DestroyVecfun(&psi_int);
	}

	for(i=H.ns; i<H.ne; i++)
		setr(&H, i, 1.0 / (1.0 + valr(&H, i) ) );
	// for plotting purposes, don't check if valr(&f, i)==0 here
	DestroyVecfun(&H);
	DestroyVecfun(&f);

}

double FormJf(Mode* ms, Geometry geo, Vec v, Vec f, double ftol, int printnewton){
	Mode m = ms[0];
	int lasing = m->lasing, Nm;
	VecGetSize(v, &Nm); Nm /= NJ(geo);
	Mat J = m->J; // for multimode, all m share same J

	if(lasing)
		ComputeGain(geo, ms, Nm); // do this before naming scratch vectors!
	// ================== name scratch vectors ================== //
	Vec vpsisq = geo->vscratch[3], // only form this later if needed
		vIpsi = geo->vscratch[2];
	Vec dfR, dfI;
	if(Nm == 1){
		dfR = geo->vscratch[0];
		dfI = geo->vscratch[1];
	}else{
		dfR = geo->vNhscratch[0];
		dfI = geo->vNhscratch[1];
	}

	// =========== linear J to compute residual ========= //
	MatRetrieveValues(J);

	int ih, jh, kh, ir, jr, jc;
	for(ih=0; ih<Nm; ih++){
		m = ms[ih];
		VecSet(dfR, 0.0);	
		VecSet(dfI, 0.0);



		LinearDerivative(m, geo, dfR, dfI, ih);
	  	SetJacobian(geo, J, dfR, -2, 0, ih);
		SetJacobian(geo, J, dfI, -2, 1, ih); 
	}


	// row derivatives already added in add placeholders!

	AssembleMat(J);
	MatMult(J, v, f);
	for(kh = 0; kh<Nm; kh++) for(ir=0; ir<2; ir++)
		VecSetValue(f, kh*NJ(geo) + Nxyzcr(geo)+ir, 0.0, INSERT_VALUES);
	// assume psi(ifix) = something + i 0. Note we do not explicitly assume what "something" is, so in principle psi can have any normalization, as long as Im psi(ifix) =0.	

	AssembleVec(f);
	double fnorm;
	VecNorm(f, NORM_2, &fnorm);

	if( printnewton ) PetscPrintf(PETSC_COMM_WORLD, "|f| = %1.6e;", fnorm);
	// no \n here to make room for timing printf statement immediately afterwards

	if(Nm==2) //DEBUG
		PetscPrintf(PETSC_COMM_WORLD, " DEBUG: |yw| c = (%g, %g)", get_c(ms[0]) / cabs(gamma_w(ms[0], geo)), get_c(ms[1]) / cabs(gamma_w(ms[1], geo)) );
		// see note under "%s at D = %g:" print statement in Newton.c

	if(fnorm < ftol )
		return fnorm;   		// TODO: deleted old integral routine. Write new one here.

	// =============== column derivatives ====================

	for(jh=0; jh<Nm; jh++){
		Mode mj = ms[jh];

		VecSet(dfR, 0.0);
		VecSet(dfI, 0.0);

		if(Nm > 1)  // hack: only recompute vpsisq if ComputeGain didn't already do it, i.e. for multimode
			VecDotMedium(geo, mj->vpsi, mj->vpsi, vpsisq, geo->vMscratch[0]);

		for(ih=0; ih<Nm; ih++){
			Mode mi = ms[ih];
			TimesI(geo, mi->vpsi, vIpsi);
			ColumnDerivative(mi, mj, geo, dfR, dfI, vIpsi, vpsisq, ih);
		}


		SetJacobian(geo, J, dfR, -1, 0, jh);
		SetJacobian(geo, J, dfI, -1, 1, jh);
	}

	//================ tensor derivatives ================

	if(lasing){
		Vec vpsibra = vpsisq; vpsisq = 0;

		for(jh=0; jh<Nm; jh++){
			Mode mj = ms[jh];

			for(jr=0; jr<2; jr++) for(jc=0; jc< geo->gN.Nc; jc++){
				VecCopy(mj->vpsi, vpsibra);
				Stamp(geo, vpsibra, jc, jr, geo->vMscratch[0]);

				VecSet(dfR, 0.0);
				for(ih=0; ih<Nm; ih++){
					Mode mi = ms[ih];
					TimesI(geo, mi->vpsi, vIpsi);
					TensorDerivative(mi, mj, geo, dfR, vpsibra, vIpsi, ih);
				}

				SetJacobian(geo, J, dfR, jc, jr, jh);
			}
		}
	}


	if(geo->interference != 0.0 && Nm == 2){
		// column interference derivative
	

		for(jh=0; jh<Nm; jh++){
			Mode mj = ms[jh];

			VecSet(dfR, 0.0);
			VecSet(dfI, 0.0);

			for(ih=0; ih<Nm; ih++){
				Mode mi = ms[ih];
				TimesI(geo, mi->vpsi, vIpsi);
				// assume vscratch[6] computed above in ComputeGain
				ColumnDerivative(mi, mj, geo, dfR, dfI, vIpsi, geo->vscratch[6], ih);
			}

			double thisc = get_c(mj), otherc = get_c(ms[1-jh]);
			AssembleVec(dfI); // hack, so VecScale can work. Usually SetJacobian does the assembling
			VecScale( dfI, otherc / thisc); // factor of two already included

			AssembleVec(dfR); 
			// hack; this is usually called in SetJacobian, but here we are not setting the dfR part. If we don't assemble here, then ColumnDerivative will complain the next time we try to insert values into dfR.
		
			SetJacobian(geo, J, dfI, -1, 1, jh);
		}

		// tensor interference derivative

		for(jc=0; jc< geo->gN.Nc; jc++){  // shifted order of loops around for convenience, but this is actually unnecessary
			Vec vpsibra = geo->vscratch[3], 
			    vcos = geo->vscratch[4], vsin = geo->vscratch[5];

			for(jh=0; jh<Nm; jh++) for(jr=0; jr<2; jr++) {	
				VecCopy( ms[1-jh]->vpsi, vcos);
				Stamp(geo, vcos, jc, jr, geo->vMscratch[0]);
				VecCopy( ms[1-jh]->vpsi, vsin);
				Stamp(geo, vsin, jc, 1-jr, geo->vMscratch[0]);



	//d/dE1R = cos E2R - sin E2I
	//d/dE1I = cos E2I + sin E2R

	//the other two are (1 <-> 2, sin <-> -sin) 


				double sign = (jh == jr ? -1.0 : 1.0);
				double costh = cos(geo->interference), sinth = sin(geo->interference);
				VecScale( vsin, sign * sinth);
				VecWAXPY( vpsibra, costh, vcos, vsin); // this function has awkward syntax, and VecAXPBY is even worse

				VecSet(dfR, 0.0);
				for(ih=0; ih<Nm; ih++){
					Mode mi = ms[ih];
					TimesI(geo, mi->vpsi, vIpsi);
					TensorDerivative(mi, ms[0], geo, dfR, vpsibra, vIpsi, ih);
				} // will have factor of c[0]^2

				double thisc = get_c(ms[0]), otherc = get_c(ms[1]);
				AssembleVec(dfR); // same hack. see above
				VecScale( dfR, otherc / thisc); // factor of two already included				

				SetJacobian(geo, J, dfR, jc, jr, jh);

			}
		}

		



		PetscPrintf(PETSC_COMM_WORLD, "DEBUG: interference code called");
	}

	AssembleMat(J);
	return fnorm;
}
