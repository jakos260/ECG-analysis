/*************************** ecgvector.c ******************************************************/
/*                                                                                    
/* ECGvector \computes the ECGvector		
/*
/* Copyright Peacs BV, Arnhem 2016-2019
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
	mexPrintf("[CIPSVECTOR,KORSVECTOR,DOWERVECTOR] = ECGvector(ECG, electrodePos,vecPos,weightsperElectrode) \n");
    mexPrintf("\n");
	mexPrintf("DESCRIPTION\n");		
    mexPrintf("ECG = for the 12 lead ECG the following signals should be used as input [aVR aVL aVf V1-V6], each row is a signal\n");
    mexPrintf("ECG = any number signals, each row is a signal\n");
    mexPrintf("electrodePos the electrode postions in the same order as the ECG is provided\n");
    mexPrintf("vecPos = the postion of the vector (usaually the mean of the ventricles, can also be time changing postion, needs one postion per sample\n");
    mexPrintf("weightsperElectrode = weight factor per electrode\n");
    mexPrintf("CIPSVECTOR the vector loop as computed from the electrode locatins\n");
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
    double *CIPSVectorLoop =NULL;
    double *KORSVectorLoop =NULL;
    double *DowerVectorLoop =NULL;
    double *PtrVectorLoop =NULL;
    
    double *elec    =NULL;
    double (*elecNorm)[3]=NULL;
    double (*vectorPos)=NULL;
    double tmp, lengthVec,x,y,z;
    
    
    double CIPSweightFixed[9] = { 3, 3, 3, 1,1,1,1,1,1}; 
//     double CIPSweightFixed[9] = { 8, 5, 4, 3,3,1,2,2,2}; 

    double *CIPSweight = NULL;
    int     xyzWeights = 0;

    
    int useVecPos = 0;
	/* Check input arguments */
	if( nrhs > 4 || nrhs < 1) 
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
         nlhs <= 3 &&
         nrhs >= 2 )
    {
        ECG =  mxGetPr(prhs[0]);
        elec = mxGetPr(prhs[1]);
        vectorPos = mxGetPr(prhs[2]);
        useVecPos = (int)mxGetM(prhs[2]) == nSamp;

        if ( nrhs == 4 &&
            ((int)mxGetM(prhs[3]) * (int)mxGetN(prhs[3])) == nSig )    
        {
//             mexPrintf("using weigths\n");
            if ( ((int)mxGetM(prhs[3]) * (int)mxGetN(prhs[3])) == nSig ) 
            {
                CIPSweight =  mxGetPr(prhs[3]); 
            }
            else if ( (int)mxGetM(prhs[3]) == nSig && 
                      (int)mxGetN(prhs[3]) == 3 ) 
            {
                CIPSweight =  mxGetPr(prhs[3]); 
                xyzWeights = 1;
            }
        }       
                        
        plhs[0]=mxCreateDoubleMatrix(nSamp,3,mxREAL);					
        CIPSVectorLoop  = mxGetPr(plhs[0]);
        
        if( nlhs > 1 &&
            nSig == 9 ) 
        {   
            plhs[1] = mxCreateDoubleMatrix(nSamp,3,mxREAL);					
            if( nlhs >= 3 )
            {
                plhs[2] = mxCreateDoubleMatrix(nSamp,3,mxREAL);		                
            }
        }
        elecNorm = mxCalloc(nElec, 3 * sizeof(double));	elecNorm || outofmem();
// mexPrintf("hallo3 %d %d %f \n",i,t,tmp);
        if ( !useVecPos ) 
        {
            for (i =  0; i < nElec; i++)
            {
                elecNorm[i][0] = elec[i            ] - vectorPos[0];
                elecNorm[i][1] = elec[i +     nElec] - vectorPos[1];		
                elecNorm[i][2] = elec[i + 2 * nElec] - vectorPos[2];
                tmp = norm(elecNorm[i]);

                if ( tmp > 0.1 )
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

                    if ( tmp > 0.1 )
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
            }
            

            for (i =  0; i < nElec; i++)
            {	

                if ( CIPSweight != NULL )
                {    
                    if ( xyzWeights )
                    {
                        tmp = ECG[i + (t * nSig)];
                        CIPSVectorLoop[t + 0 * nSamp] += elecNorm[i][0] * tmp * CIPSweight[i + 0 * nSig ];
                        CIPSVectorLoop[t + 1 * nSamp] += elecNorm[i][1] * tmp * CIPSweight[i + 1 * nSig ];
                        CIPSVectorLoop[t + 2 * nSamp] += elecNorm[i][2] * tmp * CIPSweight[i + 2 * nSig ];                        

                    }
                    else
                    {
                        tmp = ECG[i + (t * nSig)] * CIPSweight[i];
                        CIPSVectorLoop[t + 0 * nSamp] += elecNorm[i][0] * tmp;
                        CIPSVectorLoop[t + 1 * nSamp] += elecNorm[i][1] * tmp;
                        CIPSVectorLoop[t + 2 * nSamp] += elecNorm[i][2] * tmp;                        
                    }
                }
                else if ( nSig == 9 )
                {
                    tmp = ECG[i + (t * nSig)] * CIPSweightFixed[i];
                    CIPSVectorLoop[t + 0 * nSamp] += elecNorm[i][0] * tmp;
                    CIPSVectorLoop[t + 1 * nSamp] += elecNorm[i][1] * tmp;
                    CIPSVectorLoop[t + 2 * nSamp] += elecNorm[i][2] * tmp;                    
                }
                else
                {
                    tmp = ECG[i + (t * nSig)];
                    CIPSVectorLoop[t + 0 * nSamp] += elecNorm[i][0] * tmp;
                    CIPSVectorLoop[t + 1 * nSamp] += elecNorm[i][1] * tmp;
                    CIPSVectorLoop[t + 2 * nSamp] += elecNorm[i][2] * tmp;
                }
                    
// mexPrintf("hallo3 %d %d %f %f %f %f \n",i,t,tmp,elecNorm[i][0],elecNorm[i][1],elecNorm[i][2]);
            }            
        }
        mxFree(elecNorm);
        
        if( nlhs >= 2 && // kors
            nSig == 9 ) 
        {            
            KORSVectorLoop  = mxGetPr(plhs[1]);
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
                tmp = ( ECG[1 + (t * nSig)] - ECG[0 + (t * nSig)]) ;// lead I
                KORSVectorLoop[t + 0 * nSamp] += -KorsTransfer[2][6] * tmp;
                KORSVectorLoop[t + 1 * nSamp] +=  KorsTransfer[0][6] * tmp;
                KORSVectorLoop[t + 2 * nSamp] += -KorsTransfer[1][6] * tmp;
                
                tmp = ( ECG[2 + (t * nSig)] - ECG[0 + (t * nSig)]) ;// lead II
                KORSVectorLoop[t + 0 * nSamp] += -KorsTransfer[2][7] * tmp;
                KORSVectorLoop[t + 1 * nSamp] +=  KorsTransfer[0][7] * tmp;
                KORSVectorLoop[t + 2 * nSamp] += -KorsTransfer[1][7] * tmp;
            }
        }
        
        
        if( nlhs == 3 &&  // dower
            nSig == 9 ) 
        {            
            DowerVectorLoop  = mxGetPr( plhs[2] );

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
                tmp = ( ECG[1 + (t * nSig)] - ECG[0 + (t * nSig)]) ;// lead I
                DowerVectorLoop[t + 0 * nSamp] += -DowerTransfer[2][6] * tmp;
                DowerVectorLoop[t + 1 * nSamp] +=  DowerTransfer[0][6] * tmp;
                DowerVectorLoop[t + 2 * nSamp] += -DowerTransfer[1][6] * tmp;
                
                tmp = ( ECG[2 + (t * nSig)] - ECG[0 + (t * nSig)]) ;// lead II
                DowerVectorLoop[t + 0 * nSamp] += -DowerTransfer[2][7] * tmp;
                DowerVectorLoop[t + 1 * nSamp] +=  DowerTransfer[0][7] * tmp;
                DowerVectorLoop[t + 2 * nSamp] += -DowerTransfer[1][7] * tmp;
            }
        }
    }
    else
    {
        mexPrintf("nsignals: %d %d\n", nSig, nrhs);
        showHelp();
		mexErrMsgTxt("Inproper number of arguments %d ", nrhs);
    }
}

/**************************** End of ecgvector.c *****************************/
