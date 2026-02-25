#include "mex.h"
#include <math.h>


double norm(double* x);
void   vecdiff(double* x, double *y, double *z);
double det(double* x, double *y, double *z);
double solidangle(double* p, double* x1, double* x2, double* x3);

/***************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	int nv, nt, nobs;
	double *ver, *itri, *out, *obs;
	int *dim;
	int i1,i2,i3,k,n;
	double sa;
	
	/* Check for proper number of arguments. */
	if (nrhs!=3) {
		mexErrMsgTxt("3 inputs required.");
	}
	
	dim = (int*)mxGetDimensions(prhs[0]);
	nv = dim[1];
	dim = (int*)mxGetDimensions(prhs[1]);
	nt = dim[1];
	dim = (int*)mxGetDimensions(prhs[2]);
	nobs = dim[1];
	
	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix(nobs,nv, mxREAL);
	
	/* Assign pointers to each input and output. */
	ver = mxGetPr(prhs[0]);
	itri = mxGetPr(prhs[1]);
	obs = mxGetPr(prhs[2]);
	out = mxGetPr(plhs[0]);
	
	for (k=0;k<nt;k++) {
		i1 = (int)(itri[3*k]-0.5);
		i2 = (int)(itri[3*k+1]-0.5);
		i3 = (int)(itri[3*k+2]-0.5);
		
		for (n=0;n<nobs;n++) {
			sa = solidangle(obs+3*n,ver+3*i1,ver+3*i2,ver+3*i3)/3;
			out[n+nobs*i1] += sa;
			out[n+nobs*i2] += sa;
			out[n+nobs*i3] += sa;
		}
		
	}
}

/***************************************************************************/

double solidangle(double* p, double* x1, double* x2, double* x3)
/* solid angle at point p subtended by the triangle (x1,x2,x3)
   basic version (IEEE TBME 30 (1983) pp. 125-126)  */
{
	double R1[3], R2[3], R3[3];
	double lr1,lr2,lr3;
	double num, denom;
	
	vecdiff(p,x1,R1);
	vecdiff(p,x2,R2);
	vecdiff(p,x3,R3);
	
	lr1 = norm(R1);
	lr2 = norm(R2);
	lr3 = norm(R3);
	
	num = det(R1,R2,R3);
	denom=   lr1*lr2*lr3 + 
	 lr1*(R2[0]*R3[0]+R2[1]*R3[1]+R2[2]*R3[2]) +
	 lr2*(R1[0]*R3[0]+R1[1]*R3[1]+R1[2]*R3[2]) +
	 lr3*(R1[0]*R2[0]+R1[1]*R2[1]+R1[2]*R2[2]);
	
	return -2*atan2(num,denom);
}

/***************************************************************************/

void vecdiff(double* x, double *y, double *z)
{
	z[0] = y[0] - x[0];
	z[1] = y[1] - x[1];
	z[2] = y[2] - x[2];
}

double norm(double* x)
{
	return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

double det(double* x, double *y, double *z)
{
	return x[0]*(y[1]*z[2]-y[2]*z[1]) 
			 - y[0]*(x[1]*z[2]-x[2]*z[1])
			 + z[0]*(x[1]*y[2]-x[2]*y[1]);
}















