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

void TensorDerivativeCross(Mode *ms, Geometry geo, int jr, int jh, Vec df, Vec vIpsi){
// as with all cross routines, only for  Nm = 2
	// this block same as ColumnDerivativeCross
	AssembleVec(df);
	double w[2], c[2];
	dcomp yw[2];
	int i, ih;
	for(i=0; i<2; i++){
		w[i] = creal( get_w(ms[i]));
		c[i] = get_c(ms[i]);
		yw[i] = gamma_w( ms[i], geo);
	}          
	
	double *dfdpsi, *psijp1;
	VecGetArray(df, &dfdpsi);
	VecGetArray(ms[(jh+1)%2]->vpsi, &psijp1);
	Vecfun f, H;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	
	double G12 = sqr(geo->gampar) / ( sqr(geo->gampar) + sqr(w[1] - w[0]) );
	for(ih=0; ih<2; ih++){
		Complexfun psi;
		TimesI(geo, ms[ih]->vpsi, vIpsi);
		CreateComplexfun(&psi,ms[ih]->vpsi, vIpsi);

		for(i=0; i<Nxyz(geo); i++){
			if( valr(&f, i) == 0.0 ) continue;
			dcomp ksqDfHsq_ywpsi = sqr(w[ih])*geo->D * valr(&f, i) 
				* sqr(valr(&H, i) ) * yw[ih] * valc(&psi, i);
			dcomp dfdpsi_cross = ksqDfHsq_ywpsi * 2*geo->G0*G12 * c[0]*c[1]*psijp1[i+jr*Nxyz(geo)]; 
			dfdpsi[i + ih*NJ(geo)] += creal(dfdpsi_cross);
			dfdpsi[i + ih*NJ(geo) + Nxyz(geo)] += cimag(dfdpsi_cross);
		}
		DestroyComplexfun(&psi);
	}

	DestroyVecfun(&H);
	DestroyVecfun(&f);
	VecRestoreArray(df, &dfdpsi);
	VecRestoreArray(ms[(jh+1)%2]->vpsi, &psijp1);
}

void ColumnDerivativeCross(Vec dfR, Vec dfI, Vec vIpsi, Vec vhcross, Mode *ms, Geometry geo){
// for two modes near degeneracy only!
// cross term; can't put this in ColumnDerivative because need both w[2] and c[2]

	AssembleVec(dfR); AssembleVec(dfI);
	double w[2], c[2];
	dcomp yw[2];
	int i, ih;
	for(i=0; i<2; i++){
		w[i] = creal( get_w(ms[i]));
		c[i] = get_c(ms[i]);
		yw[i] = gamma_w( ms[i], geo);
	}          

	double G12 = sqr(geo->gampar) / ( sqr(geo->gampar) + sqr(w[1] - w[0]) ),
		*dfdk, *dfdc;
	VecGetArray(dfR, &dfdk);
	VecGetArray(dfI, &dfdc);
	Vecfun f, H, hcross;
	CreateVecfun(&f, geo->vf);
	CreateVecfun(&H, geo->vH);
	CreateVecfun(&hcross, vhcross);

	for(ih=0; ih<2; ih++){
		Complexfun psi;
		TimesI(geo, ms[ih]->vpsi, vIpsi);
		CreateComplexfun(&psi,ms[ih]->vpsi, vIpsi);

		for(i=0; i<Nxyzc(geo); i++){ // sequential only
			if( valr(&f, i) == 0) continue;

			dcomp ksqDfHsqhcross_ywpsi = sqr(w[ih])*geo->D * valr(&f, i) * sqr(valr(&H, i))
				* valr(&hcross, i) * yw[ih] * valc(&psi, i);
			dcomp dfdk_cross = 2.0 * ksqDfHsqhcross_ywpsi * 
					(w[ih] - w[(ih+1)%2])/sqr(geo->gampar) * G12, 
				dfdc_cross = -2.0 * ksqDfHsqhcross_ywpsi / c[ih];

			dfdk[ih*NJ(geo) + i] += creal( dfdk_cross); 
			dfdk[ih*NJ(geo) + i + Nxyzc(geo)] += cimag( dfdk_cross);
			dfdc[ih*NJ(geo) + i] += creal( dfdc_cross); 
			dfdc[ih*NJ(geo) + i + Nxyzc(geo)] += cimag( dfdc_cross);
			// same as VecSetValue in ColumnDerivative
		}
		DestroyComplexfun(&psi);
	}
	DestroyVecfun(&H);
	DestroyVecfun(&f);
	DestroyVecfun(&hcross);
	VecRestoreArray(dfR, &dfdk);
	VecRestoreArray(dfI, &dfdc);
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