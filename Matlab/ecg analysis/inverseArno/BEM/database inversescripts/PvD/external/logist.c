/****************************** logist.c **********************************/
/*                                                                         */           
/* This program calculates the logist equation for N number of points      */
/* with given R.										                              */
/* Input:  a given N and R														         */
/* Output: The logistic eqation x(n) = R*x(n-1)*(1-x(n-1))with length N    */
/*                                                                         */
/* 03-11-03   Peter van Dam                                                */
/*                                                                         */
/***************************************************************************/

#include <mex.h>
#include <math.h>

#define EPS 1e-9


/***************************************************************************/

int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}

/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	double N,R,start,(*logist); 
	int i;

  /* Check for proper number of input and output arguments. */    
  if(nrhs != 3 ) 
  		mexErrMsgTxt("3 input arguments required: 1) R, 2) start 3) N.");
  if(nlhs >2) 
   	mexErrMsgTxt("Too many output arguments.");
  /* Check data type of input argument. */

 	/* Check to make sure the last input argument is a scalar. */
	if( mxIsComplex(prhs[nrhs-1]) || mxGetN(prhs[nrhs-1])*mxGetM(prhs[nrhs-1]) != 1 ) 
  		mexErrMsgTxt("Input start must be a index number.");
	/* Get the scalar input startpnt. */
 	N = (int)mxGetScalar(prhs[nrhs-1]);
 	/* The input must be a noncomplex scalar double.*/
   if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
       !(mxGetN(prhs[0]) == 1 && mxGetM(prhs[0]) == 1)) {
     mexErrMsgTxt("Input must be a noncomplex scalar double.");
   }
   if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
       !(mxGetN(prhs[1]) == 1 && mxGetM(prhs[1]) == 1)) {
     mexErrMsgTxt("Input must be a noncomplex scalar double.");
   }
 	start = mxGetScalar(prhs[1]);
   R = mxGetScalar(prhs[0]);
  
 	
	/*Set output pointer to output matrix*/   
 	plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
 	logist = mxGetPr(plhs[0]);
 	logist[0]=start;
	for (i=1;i<N;i++)
	{
 		logist[i]=R*logist[i-1]*(1-logist[i-1]);
	}
}
