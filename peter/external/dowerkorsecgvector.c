/*************************** ecgvector.c ******************************************************/
/*                                                                                    
/* ECGvector \computes the ECGvector		
/*
/* Copyright Peacs BV, Arnhem 2016
/* SYNTAX
/*                                                                         
/***************************************************************************************************/

#include <mex.h>
#include <math.h>

/************************* Defintions and globals **************************/

// C.A. Schreurs et al. / Journal of Electrocardiology 43 (2010) 294–301
double  KorsTransfer[3][8] = {
{ -0.13,    0.05,	 0.01,	 0.14,	0.06,	0.54,	 0.38,	-0.07, },
{  0.06,   -0.02,	-0.05,	 0.06, -0.17,	0.13,	-0.07,	 0.93, },
{ -0.43,   -0.06,	-0.14,	-0.20, -0.11,	0.31,	 0.11,	-0.23, } };

 
// V1-V6 I II
 double  DowerTransfer[3][8] ={
{ -0.17,   -0.07,	 0.12,	 0.23,	0.24,	0.19,	 0.16,	-0.01, },
{  0.06,   -0.02,	-0.11,	-0.02,	0.04,	0.05,	-0.23,	 0.89, },
{ -0.23,   -0.31,	-0.25,	-0.06,	0.06,	0.11,	 0.02,	 0.1, } };


/***************************************************************************/

void showHelp(void)
{
	mexPrintf("ECG vector, computes the ECG vector loop\n\n");
    mexPrintf("All rights reserved, Copyright Peacs BV, Arnhem 2016\n");
	mexPrintf("SYNTAX\n");
	mexPrintf("[KORSVECTOR,DOWERVECTOR] = ECGvector(ECG) \n");
    mexPrintf("\n");
	mexPrintf("DESCRIPTION\n");		
    mexPrintf("ECG = for the 12 lead ECG the following signals should be used as input [aVR aVL aVf V1-V6], each row is a signal\n");
    mexPrintf("KORSVECTOR and DOWERVECTOR are computed when the 9 signals of the 12 lead ECG are used as input\n"); 
}
/***************************************************************************/
int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}

/***************************************************************************/

double norm(double r[3])
{
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  i, t , nSig=0, nSamp=0, nElec=0;
    double *ECG     =NULL;
    double *KORSVectorLoop =NULL;
    double *DowerVectorLoop =NULL;
    double *elec    =NULL;
    double (*elecNorm)[3]=NULL;
    double *pos;
    double tmp;
        
	
	/* Check input arguments */
	if( nrhs != 1  ) 
	{
		showHelp();
		mexErrMsgTxt("Inproper number of arguments ");
	}
    nSig  = (int)mxGetM(prhs[0]);
    nSamp = (int)mxGetN(prhs[0]);  	 	
    
    
    if ( nSig > 0 &&
         nlhs < 3 )
    {
        ECG =  mxGetPr(prhs[0]);
        
        if( nlhs >= 1 &&
            nSig == 9 ) 
        {      
            plhs[0]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
            KORSVectorLoop  = mxGetPr(plhs[0]);
            for (t =  0; t < nSamp; t++)
            {
                KORSVectorLoop[t  ] = 0;
                KORSVectorLoop[t +     nSamp] = 0;
                KORSVectorLoop[t + 2 * nSamp] = 0;
                // for ecg (frank etc) X = right to left, Y = head to bottom, Z = front to back
                for (i =  0; i < 6; i++)
                {	
                    tmp = ECG[i + 3 + (t * nSig)];
                    KORSVectorLoop[t + 0 * nSamp] += -KorsTransfer[2][i] * tmp;
                    KORSVectorLoop[t + 1 * nSamp] +=  KorsTransfer[0][i] * tmp;
                    KORSVectorLoop[t + 2 * nSamp] += -KorsTransfer[1][i] * tmp;
                }
                tmp = ( ECG[1 + (t * nSig)] - ECG[0 + (t * nSig)]) / 1.5;// lead I
                KORSVectorLoop[t + 0 * nSamp] += -KorsTransfer[2][6] * tmp;
                KORSVectorLoop[t + 1 * nSamp] +=  KorsTransfer[0][6] * tmp;
                KORSVectorLoop[t + 2 * nSamp] += -KorsTransfer[1][6] * tmp;
                
                tmp = ( ECG[2 + (t * nSig)] - ECG[0 + (t * nSig)]) / 1.5;// lead II
                KORSVectorLoop[t + 0 * nSamp] += -KorsTransfer[2][7] * tmp;
                KORSVectorLoop[t + 1 * nSamp] +=  KorsTransfer[0][7] * tmp;
                KORSVectorLoop[t + 2 * nSamp] += -KorsTransfer[1][7] * tmp;
            }
        }
        
        
        if( nlhs == 2 &&
            nSig == 9 ) 
        {            
            plhs[1]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
            DowerVectorLoop  = mxGetPr( plhs[1] );

            for (t =  0; t < nSamp; t++)
            {
                DowerVectorLoop[t  ] = 0;
                DowerVectorLoop[t +     nSamp] = 0;
                DowerVectorLoop[t + 2 * nSamp] = 0;
                // for ecg (frank etc) X = right to left, Y = head to bottom, Z = front to back
                for (i =  0; i < 6; i++)
                {	
                    tmp = ECG[i + 3 + (t * nSig)];
                    DowerVectorLoop[t + 0 * nSamp] += -DowerTransfer[2][i] * tmp;
                    DowerVectorLoop[t + 1 * nSamp] +=  DowerTransfer[0][i] * tmp;
                    DowerVectorLoop[t + 2 * nSamp] += -DowerTransfer[1][i] * tmp;
                }
                tmp = ( ECG[1 + (t * nSig)] - ECG[0 + (t * nSig)]) / 1.5;// lead I
                DowerVectorLoop[t + 0 * nSamp] += -DowerTransfer[2][6] * tmp;
                DowerVectorLoop[t + 1 * nSamp] +=  DowerTransfer[0][6] * tmp;
                DowerVectorLoop[t + 2 * nSamp] += -DowerTransfer[1][6] * tmp;
                
                tmp = ( ECG[2 + (t * nSig)] - ECG[0 + (t * nSig)]) / 1.5;// lead II
                DowerVectorLoop[t + 0 * nSamp] += -DowerTransfer[2][7] * tmp;
                DowerVectorLoop[t + 1 * nSamp] +=  DowerTransfer[0][7] * tmp;
                DowerVectorLoop[t + 2 * nSamp] += -DowerTransfer[1][7] * tmp;
            }
        }
    }
    else
    {
        mexErrMsgTxt("Inproper number of arguments ");
    }
}

/**************************** End of ecgvector.c *****************************/
