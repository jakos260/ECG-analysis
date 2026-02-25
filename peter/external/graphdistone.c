/*************************** graphdistone.c ******************************************************/
/*                                                                                    
/* ("GRAPHDISTONE Version 2.0 \nCalculates the shortest paths in a fully connected graph\n\n");
/* ("SYNTAX\n");
/* ("[DIST,PATHS]=graphdistone(ADJ,startpnt) \n");
/* ("DESCRIPTION\n");		
/* ("[DIST,PATH]= graphdistone(ADJ,startpnt),  Returns the distancesand (optionally) the prominent routes based on \n");
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

#define MAXNODES 5000
/************************* Defintions and globals **************************/

int		nrib,npnt;
double	maxlength,eps;
double *ADJ;
double (*DIST)		= NULL;
double (*D) 		= NULL;
char   (*fixed)		= NULL;
int    (*rib)[2]	= NULL;
double (*paths)		= NULL;


/***************************************************************************/

void showHelp(void)
{
	mexPrintf("GRAPHDISTONE \nCalculates the shortest paths in a fully connected graph\n\n");
	mexPrintf("SYNTAX\n");
	mexPrintf("[DIST,PATHS]=graphdistone(ADJ,startpnt) \n");
	mexPrintf("DESCRIPTION\n");		
	mexPrintf("[DIST,PATH]= graphdistone(ADJ,startpnt),  Returns the distancesand (optionally) the prominent routes based on \n");
	mexPrintf("                      a precomputed (ADJ) adjacency matrix of a weighted or unweighted graph. \n");
	mexPrintf("                      DIST iS the distance matrix. The distance matrix contains the distances between all nodes, while traveling along edges.\n");
	mexPrintf("                      PATH is the the prominent routes matrix. Each column (K) in PATH contains the number of times the node\n");
	mexPrintf("                      is included in all possible shortest paths starting at node K.(see van Dam PM, van Oosterom A\n");
	mexPrintf("                      Atrial Excitation Assuming Uniform Propagation. J Cardiovasc Electrophysiol. 2003;14:S166-S171)\n\n");		
	mexPrintf("See also\n- getroute or graphdist\n");
}


int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}


double epsilon(void)
 {
  volatile double e, d;
  
  for (e=1.0, d=2.0; d>1.0; e /= 2)
    d=e+1;
  return 2*e;
 }	


/***************************************************************************/
int icyc(int k,int n)
{
	if (k > n-1)
		return 0;
	else if (k <0)
		return n+k;
	else if (k>=n)
		return (k/n)-1;
	else
		return k;
}
/***************************************************************************/
static double norm(double r[3])
{
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}
/***************************************************************************/
void init_getdistmat()
{
	int i,j,k;
	int  (*trib)[2] = NULL;
	double (*tdist)	= NULL;
	
 	nrib=0;
	for ( i=0;i<npnt;i++)
	 	for (j=i;j<npnt;j++)
			if (ADJ[i*npnt+j]>0) 
				nrib++;

	rib   = mxCalloc(nrib,2*sizeof(int));   rib	  || outofmem();
	k=0;  maxlength = 0;
	/* create a list of all the ribs */
	for ( i=0;i<npnt;i++)
	 	for (j=i;j<npnt;j++)
    		if (ADJ[i*npnt+j]>0)
			{
	 		 	rib[k][0]=i;
				rib[k][1]=j;
				maxlength = maxlength+ADJ[i*npnt+j];
				k++;
			}
	fixed = mxCalloc(npnt,sizeof(char));    fixed || outofmem();
	for (i=0; i<npnt;i++)
	{	
		DIST[i]=maxlength;
		fixed[i]=0;
	}
}

/***************************************************************************/

void CalculateDists(int start)
{
	double length;
	int fixrib[MAXNODES],i,n,nfix, maxNodes=MAXNODES-1;
	init_getdistmat();
	for (i=0; i<npnt;i++)
		fixed[i]=0;
	nfix=1;
	DIST[start]=0;
	fixed[start]=1;
	/* calculate all paths between startpnt and vertices*/
	while (nfix<npnt)
	{
		length = maxlength-0.01;
		n=0;
		for (i=0;i<nrib;i++)
		{
			if ( !fixed[rib[i][1]] && fixed[rib[i][0]] &&
				  DIST[rib[i][0]]+ ADJ[rib[i][0]+rib[i][1]*npnt]<=length)
			{
				length = DIST[rib[i][0]]+ ADJ[rib[i][0]+rib[i][1]*npnt];
				fixrib[n]=i;
				n++;
			}
			else if(!fixed[rib[i][0]] && fixed[rib[i][1]] &&
						DIST[rib[i][1]]+ ADJ[rib[i][1]+rib[i][0]*npnt] <=length)
			{
				length = DIST[rib[i][1]]+ ADJ[rib[i][1]+rib[i][0]*npnt];
				fixrib[n]=i;
				n++;
			}
		}
		if (n==0)
        {
			//mexErrMsgTxt("Inconsistent graph (more than 1 part).");
            for (i = 0;i < npnt;i++)
            {
                if (!fixed[i])
                {
                    fixed[i] = 1;
                    DIST[i] = -1.0;
                }
            }
        }
		if (n>maxNodes)
			mexErrMsgTxt("Too many nodes (>1000) with the same length.");
		for (i=0;i<n;i++)
		{
			if( fixed[rib[fixrib[i]][1]] &&
				 fabs(DIST[rib[fixrib[i]][1]]+ ADJ[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < eps)      
			{
				fixed[rib[fixrib[i]][0]]=1;
				DIST[rib[fixrib[i]][0]]=length;
			}
			else if (fabs(DIST[rib[fixrib[i]][0]]+ ADJ[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < eps)      
			{
				fixed[rib[fixrib[i]][1]]=1;
				DIST[rib[fixrib[i]][1]]=length;
			}
		}
		nfix=0;
		for (i=0;i<npnt;i++)
			if (fixed[i]) 
				nfix++;
	}
	mxFree(fixed);
	mxFree(rib);

}
/***************************************************************************/

void CalculateDistsAndPaths(int start)
{
	double length;
	int fixrib[MAXNODES],i,j,n,nfix,maxNodes=MAXNODES-1,k;
	int  (*frompnt)=NULL;
	
	init_getdistmat();
	frompnt = mxCalloc(npnt,sizeof(int));    frompnt || outofmem();
	for (i=0; i<npnt;i++)
		fixed[i]=0;
	nfix=1;
	DIST[start]=0;
	fixed[start]=1;
	frompnt[start]=start;

	/* calculate all paths between startpnt and vertices*/
	while (nfix<npnt)
	{
		length = maxlength;
		n=0;
		for (i=0;i<nrib;i++)
		{
			if ( !fixed[rib[i][1]] && fixed[rib[i][0]] &&
					DIST[rib[i][0]]+ ADJ[rib[i][0]+rib[i][1]*npnt]<=length)
			{
				length = DIST[rib[i][0]] + ADJ[rib[i][0]+rib[i][1]*npnt];
				fixrib[n]=i;
				n++;
			}
			else if(!fixed[rib[i][0]] && fixed[rib[i][1]] &&
						DIST[rib[i][1]]+ ADJ[rib[i][1]+rib[i][0]*npnt] <=length)
			{
				length = DIST[rib[i][1]] + ADJ[rib[i][1]+rib[i][0]*npnt];
				fixrib[n]=i;
				n++;
			}
		}
		// if (n==0)
		// 	mexErrMsgTxt("Inconsistent graph (more than 1 part).");
		if (n>maxNodes)
			mexErrMsgTxt("Too many nodes (>1000) with the same length.");
		for (i=0;i<n;i++)
		{
			if(fixed[rib[fixrib[i]][1]] &&
				 fabs(DIST[rib[fixrib[i]][1]]+ ADJ[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < eps)      
			{
				fixed[rib[fixrib[i]][0]]=1;
				DIST[rib[fixrib[i]][0]]=length;
				frompnt[rib[fixrib[i]][0]]=rib[fixrib[i]][1];
			}
			else if (fabs(DIST[rib[fixrib[i]][0]] + ADJ[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < eps)      
			{
				fixed[rib[fixrib[i]][1]]=1;
				DIST[rib[fixrib[i]][1]]=length;
				frompnt[rib[fixrib[i]][1]]=rib[fixrib[i]][0];
			}
		}
		nfix=0;
		for (i=0;i<npnt;i++)
			if (fixed[i])
				nfix++;
	}
	for(i=0;i<npnt;i++)
	{
		j=i;
		k=0;
		paths[i+npnt*k]=i+1;	
		k=1;
		while ( j!= start )
		{
			j=frompnt[j];
			paths[i+npnt*k]=j+1;						
			k++;
		}
		if (i==start)
			paths[i+npnt*k]=j+1;
	}
	mxFree(frompnt);
	mxFree(fixed);
	mxFree(rib);

}


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  mrows,mcols,i,j,startpnt,dopath=0;
	
	npnt=0;	nrib=0;
    eps=epsilon();
	/* Check input arguments */
	if(nrhs ==1 || nrhs >2 || nlhs >2 ) 
	{
		showHelp();
		mexErrMsgTxt("Inproper number of arguments ");
		return;
	}
	else if (nrhs==2)
	{
		if (nlhs>1) dopath=1;
		/* ADJ */
	 	mrows=mxGetM(prhs[0]);		npnt = mrows;
		mcols=mxGetN(prhs[0]);  
		if( mcols !=mrows)          mexErrMsgTxt("The adjacency matrix must be square");
		if(!(mxIsDouble(prhs[0]))) 	mexErrMsgTxt("ADJ array must be of type double.");
		ADJ  = mxGetPr(prhs[0]);
		
		/* Check to make sure the last input argument is a scalar. */
		if( mxIsComplex(prhs[nrhs-1]) || mxGetN(prhs[nrhs-1])*mxGetM(prhs[nrhs-1]) != 1 ) 
			mexErrMsgTxt("Input start must be a index number.");

		startpnt = (int)mxGetScalar(prhs[nrhs-1]); startpnt--; /* Get the scalar input startpnt.*/
		plhs[0]= mxCreateDoubleMatrix(npnt,1,mxREAL);        /* Set output pointer to output matrix */
		D      = mxGetPr(plhs[0]);
        DIST=mxCalloc(npnt,sizeof(double));    DIST || outofmem();
		for (i=0; i<npnt;i++)
			DIST[i] = 0;

        if (dopath)
		{
			plhs[1] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);
			paths = mxGetPr(plhs[1]);
        }
    	if ( dopath )
    		CalculateDistsAndPaths(startpnt);
    	else 
        	CalculateDists(startpnt);
        
//    		for (i=0; i<npnt;i++)
// 			D[i]=DIST[i+startpnt*npnt];
//         
        
   		for (i=0; i<npnt;i++)
        {
            if ( DIST[i] >= maxlength )
                D[i]=-1;
            else
                D[i]=DIST[i];
        }
		mxFree(DIST);

    }
    else
        showHelp();
}

/**************************** End of graphdist.c *****************************/

