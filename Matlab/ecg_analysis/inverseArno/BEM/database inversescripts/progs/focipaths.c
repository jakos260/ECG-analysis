/*************************** graphdistone.c ******************************************************/
/*                                                                                    
/* ("fociPaths Version 2.0 \nCalculates the shortest paths in a fully connected graph\n\n");
/* ("SYNTAX\n");
/* ("[DIST]=focipaths(ADJ,foci,focitimes) \n");
/* ("DESCRIPTION\n");		
/* ("[DIST]= graphdistone(ADJ,startpnt),  Returns the distancesand (optionally) the prominent routes based on \n");
/* ("                      a precomputed (ADJ) adjacency matrix of a weighted or unweighted graph. \n");
/* ("                      DIST iS the distance matrix. The distance matrix contains the distances between all nodes, while traveling along edges.\n");
/* ("                      PATH is the the prominent routes matrix. Each column (K) in PATH contains the number of times the node\n");
/* ("                      is included in all possible shortest paths starting at node K.(see van Dam PM, van Oosterom A\n");
/* ("                      Atrial Excitation Assuming Uniform Propagation. J Cardiovasc Electrophysiol. 2003;14:S166-S171)\n\n");		
/* ("See also\n- getroute or graphdist\n");
/*
/* SEE ALSO
/*     - getroute 
/*                                                                         
/* 31-05-07   Peter van Dam                                                
/*            Medical Physics		                                          
/*                                                                         
/***************************************************************************************************/

#include <mex.h>
#include <math.h>

/************************* Defintions and globals **************************/



/***************************************************************************/

void showHelp(void)
{
	mexPrintf("focipaths \nCalculates the shortest path in a directed graph\n\n");
	mexPrintf("SYNTAX\n");
	mexPrintf("[DIST]=focipaths(ADJ,foci, focitimes) \n");
	mexPrintf("DESCRIPTION\n");		
	mexPrintf("[DIST]= focipaths(ADJ,foci, focitimes),  Returns the distances given the foci and their initial activation time based on \n");
	mexPrintf("                      a precomputed (ADJ) adjacency matrix of a (directional) weighted or unweighted graph. \n");
	mexPrintf("                      DIST iS the distance vector with the shortest distance between any of the foci and all, while traveling along edges.\n");
	mexPrintf("See also\n- getroute or graphdist\n");
}


int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}

// double min(double a, double b)
// {
//      if ( a < b)
//          return a;
//      else 
//          return b;
// }


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  tmp,i,j,k,ir;
	double *ADJ;
    int npnt;
    double *foci;
    int nFoci;
    double *fociTime;
    int nFociTime;
    int *edgeA,*edgeB,nEdge=0;
    int *fixed;
    int nfix,nprevfix;
    double *distance;
	double length,maxLength = 0;
    double minInitTime;

	npnt=0;	
	/* Check input arguments */
	if(nrhs != 3 || nlhs > 1 ) 
	{
		showHelp();
		mexErrMsgTxt("Inproper number of arguments ");
	}
	else if (nrhs == 3)
	{
		if ( mxGetM(prhs[0]) !=mxGetN(prhs[0]))          
        {
            mexErrMsgTxt("The adjacency matrix must be square");
        }
		if (!(mxIsDouble(prhs[0]))) 	
        {
            mexErrMsgTxt("ADJ array must be of type double.");
        }
		ADJ  = mxGetPr(prhs[0]);
        npnt = mxGetM(prhs[0]);
// 		if (!mxIsComplex(prhs[1]) ) 	
//         {
//             mexErrMsgTxt("foci must be of type int.");
//         }
        foci= mxGetPr(prhs[1]);
        nFoci = mxGetM(prhs[1]);
        i = mxGetN(prhs[1]);
        if ( mxGetN(prhs[1]) > nFoci )
        {
            nFoci = mxGetN(prhs[1]);
            i = mxGetM(prhs[1]);
        }
        fociTime= mxGetPr(prhs[2]);
        nFociTime = mxGetM(prhs[2]);
        j = mxGetN(prhs[2]);
        if ( nFociTime != nFoci )
        {
            nFociTime = mxGetN(prhs[2]);
            j = mxGetM(prhs[2]);
        }
        if ( nFoci == 0 )
        {
            mexErrMsgTxt("at least one focus must be given.");
        }
        if ( nFociTime != nFoci ) 
        {
            mexErrMsgTxt("number of foci and fociTimes must be equal.");
        }
        if ( i !=1 )
        {
            mexErrMsgTxt("foci must be a vector (Nx1)");
        }
        
        if ( j !=1 )
        {
            mexErrMsgTxt("fociTimes must be a vector (Nx1)");
        }
       
		plhs[0]= mxCreateDoubleMatrix(npnt,1,mxREAL);        /* Set output pointer to output matrix */
		distance = mxGetPr(plhs[0]);


        fixed = mxCalloc(npnt,sizeof(int));    fixed || outofmem();
        maxLength=0;
        for ( i=0,ir=0;i < npnt; i++,ir += npnt)
        {
            for (j=0;j<npnt;j++)
            {
                if ( ADJ[ir+j] > 0 )
                {
                    nEdge++;
                    maxLength = maxLength + ADJ[ir+j];
                }
            }
            fixed[i] = 0;
        }
       
        edgeA = mxCalloc(nEdge,sizeof(int));    edgeA|| outofmem();
        edgeB = mxCalloc(nEdge,sizeof(int));    edgeB|| outofmem();      
        k=0;
        /* create a list of all the edges*/
        for ( i=0,ir=0;i < npnt; i++,ir += npnt)
        {
            for (j=0;j<npnt;j++)
            {
                if ( ADJ[ir+j] > 0 )
                {
                    edgeA[k] = i;
                    edgeB[k] = j;
                    k++;
                }
            }
            distance[i] = maxLength;
        }

        minInitTime = maxLength;
        for (i=0;i < nFoci;i++)
        {
            foci[i] -= 1;
            if ( foci[i] >= npnt ||
                 foci[i] < 0   )
            {
                mexErrMsgTxt("foci is out of range ( > size ADJ)");
            }
            distance[(int)(foci[i])] = fociTime[i];
            if ( fociTime[i] < minInitTime )
            {
                j = i;
                minInitTime = fociTime[i];
            }
        }
       
        // only fix the initial focus. The distance of the other foci has been set but can be overruled by another focus
        nfix=0;
        for (i=0;i < nFoci;i++)
        {
            if ( fociTime[i] == minInitTime )
            {
                fixed[(int)(foci[i])] = 1;
                nfix++;
            }
        }

        nprevfix = nfix;
        /* calculate all paths between startpnt and vertices*/
        while ( nfix < npnt )
        {
            k=0;
            length = maxLength - 1.0;
            for (i=0; i < nEdge; i++ )
            {
                if ( !fixed[edgeB[i]] )
                {
                    if ( fixed[edgeA[i]] )
                    {
                        distance[edgeB[i]] = min(distance[edgeB[i]],
                                                 distance[edgeA[i]] + ADJ[edgeA[i] + edgeB[i] * npnt]);
                        length = min(length, distance[edgeB[i]] );
                    }
                    if ( k != i )
                    {
                        edgeA[k] = edgeA[i];
                        edgeB[k] = edgeB[i];
                    }
                    k++;
                }            
            }
            nEdge=k;
            nfix=0;
            for (i = 0; i < npnt;i++)
            {
                if ( distance[i] <= length )
                {
                    fixed[i] =true;
                    nfix++;
                }
            }
            if ( nfix == nprevfix )
            {
                // occurs in case on an incomplete graph
                break;
            }
            nprevfix = nfix;
        }       
		mxFree(edgeA);
        mxFree(edgeB);
        mxFree(fixed);
    }
    else
    {
        showHelp();
    }
}

/**************************** End of focipaths.c *****************************/

