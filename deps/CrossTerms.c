void ComputeHcross(Mode *ms, Geometry geo, Vec vIpsi0, Vec vIpsi1, Vec vhcross){

	VecSet(vhcross, 0.0);
	Vec vIpsi[2] = {vIpsi0, vIpsi1};
	Complexfun psi[2];
	double w[2], c[2], pval[2][2],
		G12 = sqr(geo->gampar) / ( sqr(geo->gampar) + sqr(w[1] - w[0]) );

	int i, ih;
	for(ih=0; ih<2; ih++){
		w[ih] = creal( get_w(ms[ih]));
		c[ih] = get_c(ms[ih]);
		TimesI(geo, ms[ih]->vpsi, vIpsi[ih]);
		CreateComplexfun(&psi[ih], ms[ih]->vpsi, vIpsi[ih]);
	}          
	Vecfun f, hcross;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&hcross, vhcross);
	 
	for(i=f.ns; i<f.ne; i++){
		if( valr(&f, i) == 0.0) continue;
		for(ih=0; ih<2; ih++){
			pval[ih][0] = creal( valc(&psi[ih], i) );
			pval[ih][1] = cimag( valc(&psi[ih], i) );
		}
		double val = 2*geo->G0 * G12 * c[0]*c[1]*( pval[0][0]*pval[1][0] + pval[0][1]*pval[1][1] );
		setr(&hcross, i, val);
	}

	for(ih=0; ih<2; ih++){
		DestroyComplexfun(&psi[ih]);
	}
	DestroyVecfun(&f);
	DestroyVecfun(&hcross);
}

void TensorDerivativeCross(Mode *ms, Geometry geo, int jc, int jr, int jh, Vec df, Vec vpsibra, Vec vIpsi){
// as with all cross routines, only for  Nm = 2
	// this block same as ColumnDerivativeCross
	AssembleVec(df); // switch from INSERT_VALUES to ADD_VALUES
	double w[2], c[2];
	dcomp yw[2];
	int i, ih;
	for(i=0; i<2; i++){
		w[i] = creal( get_w(ms[i]));
		c[i] = get_c(ms[i]);
		yw[i] = gamma_w( ms[i], geo);
	}          
	
	VecCopy(ms[(jh+1)%2]->vpsi, vpsibra);
	Stamp(geo, vpsibra, jc, jr, geo->vMscratch[0]);

	Vecfun f, H, psibra;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&psibra, vpsibra);
	CreateVecfun(&H, geo->vH);
	
	double G12 = sqr(geo->gampar) / ( sqr(geo->gampar) + sqr(w[1] - w[0]) );
	for(ih=0; ih<2; ih++){
		Complexfun psi;
		TimesI(geo, ms[ih]->vpsi, vIpsi);
		CreateComplexfun(&psi,ms[ih]->vpsi, vIpsi);

		for(i=f.ns; i<f.ne; i++){
			if( valr(&f, i) == 0.0 ) continue;
			dcomp ksqDfHsq_ywpsi = sqr(w[ih])*geo->D * valr(&f, i) 
				* sqr(valr(&H, i) ) * yw[ih] * valc(&psi, i);
			dcomp dfdpsi_cross = ksqDfHsq_ywpsi * 2*geo->G0*G12 * c[0]*c[1]*valr(&psibra, i); 

			int ir = i/Nxyzc(geo);
			VecSetValue(df, i + ih*NJ(geo), ir? cimag(dfdpsi_cross) : creal(dfdpsi_cross), ADD_VALUES );
		}
		DestroyComplexfun(&psi);
	}
	DestroyVecfun(&H);
	DestroyVecfun(&f);
}

void ColumnDerivativeCross(Vec dfR, Vec dfI, Vec vIpsi, Vec vhcross, Mode *ms, Geometry geo){
// for two modes near degeneracy only!
// cross term; can't put this in ColumnDerivative because need both w[2] and c[2]

	AssembleVec(dfR); AssembleVec(dfI); // switch from INSERT_VALUES to ADD_VALUES
	double w[2], c[2];
	dcomp yw[2];
	int i, ih;
	for(i=0; i<2; i++){
		w[i] = creal( get_w(ms[i]));
		c[i] = get_c(ms[i]);
		yw[i] = gamma_w( ms[i], geo);
	}          

	double G12 = sqr(geo->gampar) / ( sqr(geo->gampar) + sqr(w[1] - w[0]) );
	Vecfun f, H, hcross;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&hcross, vhcross);

	for(ih=0; ih<2; ih++){
		Complexfun psi;
		TimesI(geo, ms[ih]->vpsi, vIpsi);
		CreateComplexfun(&psi,ms[ih]->vpsi, vIpsi);

		for(i=f.ns; i<f.ne; i++){
			if( valr(&f, i) == 0) continue;

			dcomp ksqDfHsqhcross_ywpsi = sqr(w[ih])*geo->D * valr(&f, i) * sqr(valr(&H, i))
				* valr(&hcross, i) * yw[ih] * valc(&psi, i);
			dcomp dfdk_cross = 2.0 * ksqDfHsqhcross_ywpsi * 
					(w[ih] - w[(ih+1)%2])/sqr(geo->gampar) * G12, 
				dfdc_cross = -2.0 * ksqDfHsqhcross_ywpsi / c[ih];

			int ir = i / Nxyzc(geo);
			VecSetValue(dfR, ih*NJ(geo) + i, ir? cimag(dfdk_cross) : creal(dfdk_cross), ADD_VALUES);
			VecSetValue(dfI, ih*NJ(geo) + i, ir? cimag(dfdc_cross) : creal(dfdc_cross), ADD_VALUES);
		}
		DestroyComplexfun(&psi);
	}
	DestroyVecfun(&H);
	DestroyVecfun(&f);
	DestroyVecfun(&hcross);
	AssembleVec(dfR); AssembleVec(dfI);
}

void AddCrossTerm(Geometry geo, Vec vhcross){
// need to do this separated from ComputeGain because ComputeGain done before naming scratch vectors
	Vecfun f, H, hcross;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&hcross, vhcross);

	int i;
	for(i=H.ns; i<H.ne; i++){
		if( valr(&f, i) == 0.0) continue;
		double val = 1.0/valr(&H, i) + valr(&hcross, i);
		setr( &H, i, 1.0/val);
	}

	DestroyVecfun(&f);
	DestroyVecfun(&H);
	DestroyVecfun(&hcross);
}