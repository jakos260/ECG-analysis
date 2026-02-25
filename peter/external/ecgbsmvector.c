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
	mexPrintf("[CIPSVECTOR] = ECGvector(ECG, electrodePos,vecPos) \n");
    mexPrintf("\n");
	mexPrintf("DESCRIPTION\n");		
    mexPrintf("ECG = any ECG, for the 12 lead ECG the following signals should be used as input [aVR aVL aVf V1-V6], each row is a signal\n");
    mexPrintf("ECG = any number signals, each row is a signal\n");
    mexPrintf("electrodePos the electrode postions in the same order as the ECG is provided\n");
    mexPrintf("vecPos = the postion of the vector (usaually the mean of the ventricles, can also be time changing postion, needs one postion per sample\n");
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
    int useVecPos = 0;
    double *ECG     =NULL;
    double *CIPSVectorLoop =NULL;
    
    double *elec    =NULL;
    double (*elecNorm)[3]=NULL;
    double *elecNormLen=NULL;

    double (*vectorPos)=NULL;
    double tmp, lengthVec,x,y,z;
    double maxLen = 0;
    
	/* Check input arguments */
	if( nrhs != 3) 
	{
        mexPrintf("nsignals hier: %d %d\n", nSig, nrhs);
        showHelp();
        mexErrMsgTxt("Inproper number of arguments ");
	}

    
    nSig  = (int)mxGetM(prhs[0]);
    nSamp = (int)mxGetN(prhs[0]);  	 	
    nElec = (int)mxGetM(prhs[1]);
    
   
    if ( nSig > 0 &&
         nSig == nElec && 
         (int)mxGetN(prhs[1]) == 3  &&
        ((int)mxGetM(prhs[2]) == 1 || (int)mxGetM(prhs[2]) == nSamp ) &&
         (int)mxGetN(prhs[2]) == 3 && 
         nlhs <= 1 &&
         nrhs >= 2 )
    {
        ECG =  mxGetPr(prhs[0]);
        elec = mxGetPr(prhs[1]);
        vectorPos = mxGetPr(prhs[2]);
        useVecPos = (int)mxGetM(prhs[2]) == nSamp;

                      
        plhs[0]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
        CIPSVectorLoop  = mxGetPr(plhs[0]);
                
        elecNorm = mxCalloc(nElec, 3 * sizeof(double));	elecNorm || outofmem();
        elecNormLen = mxCalloc(nElec, sizeof(double));	elecNormLen || outofmem();
// mexPrintf("hallo3 %d %d %f \n",i,t,tmp);
        
        for (i =  0; i < nElec; i++)
        {
            elecNorm[i][0] = elec[i            ] - vectorPos[0];
            elecNorm[i][1] = elec[i +     nElec] - vectorPos[1];		
            elecNorm[i][2] = elec[i + 2 * nElec] - vectorPos[2];
            tmp = norm(elecNorm[i]);
            elecNormLen[i] = tmp;
//             if ( tmp > 0.1 )
//             {
//                 elecNorm[i][0] /= tmp;
//                 elecNorm[i][1] /= tmp;
//                 elecNorm[i][2] /= tmp;
//             }
//             else
//             {
//                 elecNorm[i][0] = 0;
//                 elecNorm[i][1] = 0;
//                 elecNorm[i][2] = 0;
//             }  
        }
        

        for (t =  0; t < nSamp; t++)
        {
            CIPSVectorLoop[t  ] = 0;
            CIPSVectorLoop[t +     nSamp] = 0;
            CIPSVectorLoop[t + 2 * nSamp] = 0;

            if ( useVecPos ) 
            {
                for (i =  0; i < nElec; i++)
                {
                    elecNorm[i][0] = elec[i            ] - vectorPos[t            ];
                    elecNorm[i][1] = elec[i +     nElec] - vectorPos[t + nSamp    ];		
                    elecNorm[i][2] = elec[i + 2 * nElec] - vectorPos[t + nSamp * 2];
                    tmp = norm(elecNorm[i]);
                    elecNormLen[i] = tmp;

//                     if ( tmp > 0.1 )
//                     {
//                         elecNorm[i][0] /= tmp;
//                         elecNorm[i][1] /= tmp;
//                         elecNorm[i][2] /= tmp;
//                     }
//                     else
//                     {                                                
//                         elecNorm[i][0] = 0;
//                         elecNorm[i][1] = 0;
//                         elecNorm[i][2] = 0;
//                     } 
                }
            }
            

            for (i =  0; i < nElec; i++)
            {	
                tmp = ECG[i + (t * nSig)] * elecNormLen[i];
                    
                CIPSVectorLoop[t + 0 * nSamp] += (elecNorm[i][0] * tmp);
                CIPSVectorLoop[t + 1 * nSamp] += (elecNorm[i][1] * tmp);
                CIPSVectorLoop[t + 2 * nSamp] += (elecNorm[i][2] * tmp);
            }     
            tmp = sqrt( CIPSVectorLoop[t + 0 * nSamp] * CIPSVectorLoop[t + 0 * nSamp] + 
                        CIPSVectorLoop[t + 1 * nSamp] * CIPSVectorLoop[t + 1 * nSamp] +
                        CIPSVectorLoop[t + 2 * nSamp] * CIPSVectorLoop[t + 2 * nSamp] );   
            if ( tmp > maxLen )
            {
                maxLen = tmp;
            }
        }
        for (t =  0; t < nSamp; t++)
        {
            CIPSVectorLoop[t + 0 * nSamp] /= maxLen;
            CIPSVectorLoop[t + 1 * nSamp] /= maxLen;
            CIPSVectorLoop[t + 2 * nSamp] /= maxLen;
        }
        mxFree(elecNorm);
    }
    else
    {
        mexPrintf("nsignals: signals %d inputs %d elec %d\n", nSig, nrhs,nElec);
        showHelp();
		mexErrMsgTxt("Inproper number of arguments");
    }
}


// mexPrintf("hallo3 %d %d %f %f %f %f \n",i,t,tmp,elecNorm[i][0],elecNorm[i][1],elecNorm[i][2]);

/**************************** End of ecgvector.c *****************************/
