#include "headers.h"

dcomp pmlval(int i, int* N, int* Npml, double* h, int LowerPML, int k){
// k = 1 for muinv, k = 0 for eps
	Point p(i, Grid(N, 3, 2) );

	double omega = 1.0, d[3], val[2] = {0.0, 0.0};			
	for(int j=0; j<3; j++) 
	if( Npml[j]==0) d[j] = 0.0;
	else{

int 	lower = Npml[j]-p.x(j),
upper = p.x(j) - (N[j]-Npml[j]-1-(k? j!=p.c(): j==p.c()));

d[j] = LowerPML*(lower>0)*(lower - 0.5*(k==0? j==p.c() : j!=p.c()))/Npml[j]
+ (upper>0)*(upper - 0.5*(k==0? j==p.c() : j!=p.c()))/Npml[j];
	}
	

	double sigma[3];
	for(int j=0; j<3; j++) 
		sigma[j] = Npml[j] == 0 ? 0 : -3.0/4*log(1e-25)/(Npml[j]*h[j]);



	for(int j=0; j<3; j++) 
		val[0] += sigma[j]* sigma[(j+1)%3]*
			     sqr(d[j]*d[(j+1)%3]) * (j==(p.c()+1)%3?-1 : 1);

	val[0] = val[0]/sqr(omega)+1;

	for(int j=0; j<3; j++) 
		val[1] += sigma[j] *sqr(d[j])*(j==p.c()? -1:1);

	val[1] += sigma[0]* sigma[1]*sigma[2]
			*sqr(d[0]*d[1]*d[2])/sqr(omega);
	val[1] /= omega * (k==0? 1: -1);

	for(int l=0; l<2; l++) for(int j=(k? 1:0); j<(k?3:1); j++) 
		val[l] /= 1 + sqr( sigma[(j+p.c())%3] *sqr(d[(p.c()+j)%3])/omega);

	return val[0]+ComplexI*val[1];


}
