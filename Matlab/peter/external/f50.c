
/*	f50.c: 50 Hz notch filter	*/

/* source (?):

        Mortara D.W.     Digital filters for ecg signals.
	Proc. Computers in Cardiology, 511-513.
*/
		

#include <string.h>
#include "math.h"
#include "mex.h"

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	double	*s,*S;
	double	meandc=0, *mean=NULL, sampfreq=1000,nmean;
	double	alfa=0.9, beta, gamma, corr, tmp;
	int	i,k, imean;
	int	irow, icol, nrow=0, ncol=0;

	/* Check input arguments */
	if(nrhs >3 ||nrhs==0) 
	{
	    mexPrintf("sharpness of the filter (alpha: 0.5 - 1; 1 is sharpest): ");
		mexErrMsgTxt("Maximally 3 input parameters are required: 1) signal, 2) sample frequency 3) alpha");
	}
	if(nrhs >1)
	{
		if( mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1 )
			mexErrMsgTxt("second parameter should be a scalar (sample frequency [Hz]) ");
		else
			sampfreq=mxGetScalar(prhs[1]);

		if(nrhs ==3)
		{
			if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )
				mexErrMsgTxt("Third parameter should be a scalar: sharpness of the filter (alpha: 0.5 - 1; 1 is sharpest): ");
			else
				alfa=mxGetScalar(prhs[2]);
		}		
	}
	S=mxGetPr(prhs[0]);
	nrow=mxGetM(prhs[0]);
	ncol=mxGetN(prhs[0]); 
	
	mexPrintf("sampleFreqency %6.1f     alpha %5.3f \n",sampfreq,alfa);						
	plhs[0]=mxCreateDoubleMatrix(nrow,ncol,mxREAL);
	s= mxGetPr(plhs[0]);		
	for (i=0;i<nrow*ncol;i++)	
		s[i]=S[i];
	
	/*	determin coefficients	*/
	nmean=(20*sampfreq/1000);	/* number of samples in 20 ms */
	gamma=1-(1-alfa)*(1-alfa)/(nmean*nmean*(alfa-2*alfa*alfa));
	beta=((gamma+1)-sqrt(gamma-1)*sqrt(gamma+3))/2;
	corr=(beta/alfa)*(1+alfa)/(1+beta);
	mean=mxCalloc(ceil(nmean),sizeof(double));
	if (mean==NULL)
		mexErrMsgTxt("Cannot allocate mean\n");
	
	/*	do the filtering	*/
	
	for (irow=0; irow<nrow; irow++)
	{
		meandc=0;
		for (k=0; k<nmean; k++)
        {
			mean[k]=0;
        }
        imean = 0;
		k = (10*nmean/(1-alfa));
		if (k>ncol-1) k=ncol-1;
		for (icol=k; icol>=0; icol--)	/* initiate imean in reverse order */
		{
            imean = ((icol+1) - (icol / nmean) * nmean) - 1;			
			mean[imean]=(1-alfa)*s[irow+icol*nrow]+alfa*mean[imean];
			meandc=(1-beta)*s[irow+icol*nrow]+beta*meandc;
		}
		imean=0;
		for (icol=0; icol<ncol; icol++)	/* the actual filtering		   */
		{			
			mean[imean]=(1-alfa)*s[irow+icol*nrow]+alfa*mean[imean];
			meandc=(1-beta)*s[irow+icol*nrow]+beta*meandc;
			s[irow+icol*nrow]=corr*(s[irow+icol*nrow]-mean[imean])+meandc;
            imean++;
			if (imean>=nmean)
				imean=0;
		}
	}
	mxFree(mean);
}  
