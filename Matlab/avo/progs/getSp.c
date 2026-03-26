/*************************** getsp.c    ***********************************/
/*                                                                         */           
/* This program calculates the shortest paths between one vertex and all   */
/* the other vertices of a triangulated object.                            */
/* The program creates matrix of 1 by number of vertices.				   */
/* Input:  An adjaceny matrix (file) with the distances between the        */
/*         adjacent vertices                                               */
/* Output: A matrix with all shortest paths between the all vertices.      */
/*                                                                         */
/* 20-2-14   Peter van Dam                                                 */
/*			  Peacs						                                   */  
/*                                                                         */
/***************************************************************************/

#include <mex.h>
#include <math.h>

/************************* Definitions and global **************************/

int		mNpnt		= 0;
int 	mMaxt		= 0;
int 	mMode		= 1;
int 	mSourceVersion  = 0;
int    *mNeighDep	= NULL;
double *mS  		= NULL;
const double *mADJ		= NULL;
const double *mDep		= NULL;
const double *mRep		= NULL;
const double *mRepParams	= NULL;


/***************************************************************************/
void showHelp(void)
{
    mexPrintf("getSp \n compute the source strength\n");
    mexPrintf("SYNTAX\n");
    mexPrintf("S = getSp(dep, rep, maxt,ADJ,mode,repparams) \n");

    mexPrintf("DESCRIPTION\n");
    mexPrintf("S = getSp(dep, rep, maxt,ADJ,mode,repparams),  \n");
	mexPrintf("mode: 1) dep only \n");
	mexPrintf("		 2) S according to ECGsim, requires plateau slope and rep slope\n");
	mexPrintf("		 3) S using the cumsum (avo), requires a) initial repslope, b) plateau slope and c)rep slope\n");
	mexPrintf("maxt: maximum time, is the number of columns of S \n");
    mexPrintf("dep : the depolarization time \n");
	mexPrintf("ADJ : adjacency matrix with graphdist(ITRI) \n");
	mexPrintf("rep : the depolarization time \n");
	mexPrintf("repparams: This depends on the mode \n");
	mexPrintf("		  mode == 2	 plateau slope and rep slope \n");
	mexPrintf("		  mode == 3	 initial repslope, b) plateau slope and c)rep slope \n");
}


int outofmem(void)
{
    mexErrMsgTxt("Out of memory.");
    return 1;
}

void computedSrowDep(int ir,double * rowS )
{
	int i, ti, curDepI, maxToZeroT,mint,maxt;
	int Nneigh;
	double t, slope,factor; 
	const double *rowAdj;
    
	rowAdj 		= &mADJ[ir * mNpnt];
	
	curDepI = floor( mDep[ir] + 0.499 );
	maxToZeroT 	= curDepI;

	Nneigh = 0;
	for (i=0; i < mNpnt; i++ )
	{			
		if ( rowAdj[i] > 0 )
		{
			mNeighDep[Nneigh] = floor( mDep[i] + 0.4999);

			if ( mNeighDep[Nneigh] > maxToZeroT )
			{
				maxToZeroT = mNeighDep[Nneigh];
			}
			Nneigh++;
		}		
	}
	if ( maxToZeroT >= mMaxt )
	{
		maxToZeroT = mMaxt - 1;
	}
	
    for (ti = 0; ti <= maxToZeroT; ti++ )
    {
        rowS[ti] = 0;
    }
	for ( ; ti < mMaxt; ti++ )
	{
		rowS[ti] = 1.0;
	}

    factor = 1.0 / (double) Nneigh;
    for (i = 0; i < Nneigh; i++ )
    {
		mint = ( curDepI < mNeighDep[i] ) ? curDepI : mNeighDep[i];
		if ( mint < mMaxt )
		{
			maxt = ( curDepI > mNeighDep[i] ) ? curDepI : mNeighDep[i];
			if ( maxt > mMaxt )
			{
				maxt = mMaxt;
			}

			slope = factor / (double)(maxt - mint);
			t = 0;        
			for ( ti = mint; ti < maxt; ti++, t += 1.0 )
			{
				rowS[ti] += slope * t;            
			}

			for ( ; ti <= maxToZeroT; ti++)
			{
				rowS[ti] += factor;
			}
		}
    }
}



void getS( )
{     
	int t, ir;
	double tdep, trep, maxAmpl, sumY, vali;
	double *rowS;
	double *y;

	y = mxCalloc(mMaxt, sizeof(double));	y || outofmem();
	
	mNeighDep = mxCalloc(mNpnt, sizeof(int));	mNeighDep || outofmem();
	
    for ( ir = 0; ir < mNpnt; ir++ )
    {		
        rowS = &mS[ ir * mMaxt ];
		
        tdep    = -mDep[ir];
        trep    = -mRep[ir];
        maxAmpl = -1e6;

        if ( 0 ) // will create a spiky signal due to the steep upslope
        {
            for ( t = 0; t < mMaxt; t++, tdep += 1.0, trep += 1.0 )
            {
                rowS[t] = 1.0 /( 1.0 +  exp(-2.0  * tdep ) );
            }
            tdep    = -mDep[ir];
        }
        else
        {		
            computedSrowDep( ir, rowS );
        }

        if ( mMode > 1 )
        {
            if ( mSourceVersion == 1 )
            {
                for ( t = 0; t < mMaxt; t++, tdep += 1.0, trep += 1.0 )
                {
                    vali =( 1.0 / ( 1.0 + exp( mRepParams[0] * trep ) ) ) *
                          ( 1.0 / ( 1.0 + exp( mRepParams[1]     * trep ) ) );

                    rowS[t] = vali * rowS[t];

                    vali *= 1.0 /( 1.0 +  exp(-2.0 * tdep ) );
                    maxAmpl = maxAmpl > vali ? maxAmpl : vali;
                }
            }
            else // mSourceVersion == 0
            {
                sumY = 0;

                tdep    = -mDep[ir];
                trep    = -mRep[ir];

                for ( t = 0; t < mMaxt; t++, tdep += 1.0, trep += 1.0 )
                {
                    sumY +=   mRepParams[0] +
                                ( 1 / ( 1 + exp( mRepParams[1] * trep ) ) ) *
                                ( 1 / ( 1 + exp( mRepParams[2]     * trep ) ) );

                    y[t] = sumY; // 1 - cumsum
                }
                for (t=0; t < mMaxt; t++ )
                {
                    rowS[t] *= 1- (y[t]/sumY);      // scale between 0 - 1
                                                    // multiply dep and rep S = S * (1-Y)
                    maxAmpl = maxAmpl > rowS[t] ? maxAmpl : rowS[t];
                }


            }
            // For the future when dealing with ischemic areas
            //        maxAmpl *= ampl.at(ir);// devide by maxAmpl to normalize
            for (t=0; t < mMaxt; t++ ) // scale between 0 - 1
            {
                rowS[t] /= maxAmpl;
            }
        }
    }

	mxFree( y );
	mxFree( mNeighDep );
}

/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{		
    /* Check input arguments */
    if (nrhs < 4)
    {
        showHelp();
		mexErrMsgTxt("at least 4 input arguments required (mode ==1)");
        return;
    }
    if( nrhs > 6 )
    {
        showHelp();
        mexErrMsgTxt("maximum 6 input arguments required");
    }
    if(nlhs > 1)
    {
        showHelp();
        mexErrMsgTxt("only one output arguments is defined. ");
    }
	
	// mode
	if( mxGetM(prhs[0]) * mxGetN(prhs[0]) != 1 )
	{
		showHelp();
		mexErrMsgTxt("mode is a single number 1: dep only, 2,3: dep and rep");
	}
	mMode = (int)mxGetScalar(prhs[0]);	
	
	// maxt
	if( mxGetM(prhs[1]) > 1  || mxGetN(prhs[1]) > 1 )
	{
		showHelp();
		mexErrMsgTxt("maxt is a single number");
	}
	mMaxt = (int)mxGetScalar(prhs[1]);
	
// ADJ
		if ( mxGetN(prhs[2]) == mxGetM(prhs[2]) ) //ADJ
		{
			mADJ = mxGetPr(prhs[2]);
			mNpnt = (int)mxGetM(prhs[2]);
		}
        else
        {
            showHelp();
            mexErrMsgTxt("size ADJ = [n, n]");
        }

	// dep
	if( !(mxIsDouble(prhs[3])) && 
		 (int)(mxGetM(prhs[3]) * mxGetN(prhs[3])) != mNpnt )
	{
		showHelp();
		mexErrMsgTxt("dep(n,1) of type double");
	}
	mDep  = mxGetPr(prhs[3]);
				
	if ( mMode > 1 && nrhs != 6 )
	{
		showHelp();
		mexErrMsgTxt("number of parameters does not match with the selected mode");
	}

	// rep
	if( (!mxIsDouble(prhs[4]) && (int)(mxGetM(prhs[4]) * mxGetN(prhs[4])) != mNpnt) )
	{
		showHelp();
		mexErrMsgTxt("rep(n,1) of type double when mode not == 1");
	}
	mRep  = mxGetPr(prhs[4]);
			
        
		
	if ( mMode == 1 )
	{
		// no rep params required
	}
	else if ( mMode == 2 )
	{
		mSourceVersion = 1; // ECGsim
		if ( mxGetN(prhs[5]) * mxGetM(prhs[5]) == 2 &&
			 mxIsDouble(prhs[5]) ) 
		{
			mRepParams = mxGetPr(prhs[5]);		
		}			
		else
		{
			showHelp();
			mexErrMsgTxt("not the right repolarization parameters");
		}
	}
	else if ( mMode == 3 )
	{
		mSourceVersion = 0; // AVO
		if ( mxGetN(prhs[5]) * mxGetM(prhs[5]) == 3 &&
			 mxIsDouble(prhs[5]) ) 
		{
			mRepParams = mxGetPr(prhs[5]);		
		}			
		else
		{
			showHelp();
			mexErrMsgTxt("not the right repolarization parameters");
		}
	}
	else
	{
		showHelp();
		mexErrMsgTxt("no repolarization parameters or the wrong mode");
	}					
		        
	plhs[0] = mxCreateDoubleMatrix(mMaxt,mNpnt,mxREAL);					
	mS = mxGetPr( plhs[0] );
	getS();
}

/**************************** End of distgraph.c *****************************/


