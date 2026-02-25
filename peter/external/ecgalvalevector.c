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
    mexPrintf("All rights reserved, Copyright Peacs BV, Arnhem 2019\n");
	mexPrintf("SYNTAX\n");
	mexPrintf("[VCG,MEANTSI] = ECGvector(ECG, electrodePos,startPos) \n");
    mexPrintf("\n");
	mexPrintf("DESCRIPTION\n");		
    mexPrintf("ECG = for the 12 lead ECG the following signals should be used as input [aVR aVL aVf V1-V6], each row is a signal\n");
    mexPrintf("ECG = any number signals, each row is a signal\n");
    mexPrintf("electrodePos the electrode postions in the same order as the ECG is provided\n");
    mexPrintf("startPos = the postion of the meanTSI (usaually the mean of the ventricles, can also be time changing postion, needs one postion per sample\n");
    mexPrintf("VCG the vector loop as computed from the electrode locations\n");
    mexPrintf("MEANTSI the position of the meanTSI over time from the VCG\n"); 
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
    
    double *VCG =NULL;    
    double *meanTSI =NULL;
    
    double *elec    =NULL;
    double (*elecNorm)[3]=NULL;
    double (*startPos)=NULL;
    double direction[3];
    double tmp, maxVcg;
    double rmsECG,maxRmsECG;
    
    double minCosAngle = 0.5;
    
   
    
	/* Check input arguments */
	if( nrhs != 3 &&
        nlhs != 2 ) 
	{
        mexPrintf("nsignals hier: %d %d\n", nSig, nrhs);
        showHelp();
        mexErrMsgTxt("Inproper number of input arguments (3 input required, 2 output)");
	}

    
    nSig  = (int)mxGetM(prhs[0]);
    nSamp = (int)mxGetN(prhs[0]);  	 	
    nElec = (int)mxGetM(prhs[1]);
    
   
    if ( nSig > 0 &&
         nSig == nElec && 
         (int)mxGetN(prhs[1]) == 3  &&
        ((int)mxGetM(prhs[2]) == 1 || (int)mxGetM(prhs[2]) == nSamp ) &&
         (int)mxGetN(prhs[2]) == 3 && 
         nlhs <= 3 &&
         nrhs >= 2 )
    {
        ECG =  mxGetPr(prhs[0]);
        elec = mxGetPr(prhs[1]);
        startPos = mxGetPr( prhs[2] );        
                        
        plhs[0]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
        VCG  = mxGetPr(plhs[0]);
        plhs[1]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
        meanTSI  = mxGetPr(plhs[1]);
                
        elecNorm = mxCalloc(nElec, 3 * sizeof(double));	elecNorm || outofmem();
        

        t=0;
        VCG[t            ] = 0;
        VCG[t +     nSamp] = 0;
        VCG[t + 2 * nSamp] = 0;
        meanTSI[t            ] = startPos[0];
        meanTSI[t +     nSamp] = startPos[1];
        meanTSI[t + 2 * nSamp] = startPos[2];
        
        
        maxRmsECG = 0;
        maxVcg = 0;
        for (t =  1; t < nSamp; t++)
        {
            VCG[t            ] = 0;
            VCG[t +     nSamp] = 0;
            VCG[t + 2 * nSamp] = 0;

            // Vr, Vl and Vf are the first 3 electrodes
            for (i =  0; i < 3; i++)
            {
                elecNorm[i][0] = 0;
                elecNorm[i][1] = (elec[i +     nElec] - meanTSI[(t-1) + nSamp    ]);		
                elecNorm[i][2] = (elec[i + 2 * nElec] - meanTSI[(t-1) + nSamp * 2]);                                             
                tmp = norm( elecNorm[i] ) / 2.0;
//                 elecNorm[i][0] /= tmp;
                elecNorm[i][1] /= tmp;
                elecNorm[i][2] /= tmp;
            }
            for (i =  3; i < nElec; i++)
            {
                elecNorm[i][0] = elec[i            ] - meanTSI[(t-1)            ];
                elecNorm[i][1] = elec[i +     nElec] - meanTSI[(t-1) + nSamp    ];		
                elecNorm[i][2] = elec[i + 2 * nElec] - meanTSI[(t-1) + nSamp * 2];  
                tmp = norm( elecNorm[i] );
                if ( tmp > 1e-9 )
                {
                    elecNorm[i][0] /= tmp;
                    elecNorm[i][1] /= tmp;
                    elecNorm[i][2] /= tmp;
                }
                else
                {
                    elecNorm[i][0] = 0;
                    elecNorm[i][1] = 0;
                    elecNorm[i][2] = 0;
                }
            }               
            
            rmsECG = 0; 
            for (i = 0; i < nElec; i++)
            {	
                tmp = ECG[i + (t * nSig)];                
                rmsECG += (tmp * tmp);                                    
                              
                VCG[t + 0 * nSamp] += (elecNorm[i][0] * tmp);                
                VCG[t + 1 * nSamp] += (elecNorm[i][1] * tmp);
                VCG[t + 2 * nSamp] += (elecNorm[i][2] * tmp);
            }   
//             mexPrintf("elecW %f %f %f \n",elecW[0],elecW[1],elecW[2]);

            rmsECG = sqrt(rmsECG) / (double)nElec;
            if ( rmsECG > maxRmsECG )
            {
                maxRmsECG = rmsECG;
            }
            tmp =     sqrt( (VCG[t + 0 * nSamp] * VCG[t + 0 * nSamp]) + 
                            (VCG[t + 1 * nSamp] * VCG[t + 1 * nSamp]) + 
                            (VCG[t + 2 * nSamp] * VCG[t + 2 * nSamp]) );
            if ( tmp > maxVcg ) 
            {
                maxVcg = tmp;
            }

            direction[0] = VCG[t + 0 * nSamp];
            direction[1] = VCG[t + 1 * nSamp];
            direction[2] = VCG[t + 2 * nSamp];
            tmp = norm(direction);

            meanTSI[t            ] = meanTSI[(t-1)            ] + 0.7 * direction[0] / tmp;
            meanTSI[t +     nSamp] = meanTSI[(t-1) +     nSamp] + 0.7 * direction[1] / tmp;
            meanTSI[t + 2 * nSamp] = meanTSI[(t-1) + 2 * nSamp] + 0.7 * direction[2] / tmp;
        }
        maxRmsECG = 1.0;
        for (t =  1; t < nSamp; t++)
        {
            VCG[t + 0 * nSamp] *= (maxRmsECG / maxVcg); 
            VCG[t + 1 * nSamp] *= (maxRmsECG / maxVcg); 
            VCG[t + 2 * nSamp] *= (maxRmsECG / maxVcg);  
        }
        mxFree(elecNorm);                
    }
    else
    {
        mexPrintf("nsignals: %d %d\n", nSig, nrhs);
        showHelp();
		mexErrMsgTxt("Inproper number of arguments");
    }
}

// mexPrintf("hallo3 %d %d %f %f %f %f \n",i,t,tmp,elecNorm[i][0],elecNorm[i][1],elecNorm[i][2]);


/**************************** End of ecgvector.c *****************************/
