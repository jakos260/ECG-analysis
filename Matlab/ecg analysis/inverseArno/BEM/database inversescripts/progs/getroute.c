/*************************** getroute.c ***********************************/
/*                                                                         */           
/* This program calculates the shortest paths between one vertex and all   */
/* the other vertices of a triangulated object.                            */
/* The program creates matrix of 1 by number of vertices.						*/
/* Input:  An adjaceny matrix (file) with the distances between the        */
/*         adjacent vertices                                               */
/* Output: A matrix with all shortest paths between the all vertices.      */
/*                                                                         */
/* 10-11-04   Peter van Dam                                                */
/*            Medical Physics		                                          */
/*			     Vitatron						                                    */  
/*                                                                         */
/***************************************************************************/

#include <mex.h>
#include <math.h>

#define EPS 1e-9
#define MAXDHKVANRIB 12
#define MAXNODES 1000
/************************* Defintions and globals **************************/

int		nrib,npnt,ndhk,t;
double	maxlength;
double	length;
int	 (*dhk)[3]		=NULL;
double (*ver)			= NULL;
double (*D)				= NULL;
char	 (*fixed)		= NULL;
int	 (*rib)[2]		= NULL;
int	  startPnt,endPnt;
	int  (*frompnt)	=NULL;
/***************************************************************************/	
	
void showHelp()
{
	mexPrintf("Getroute calculates the shortest path between a start and end node (index).\n\n");	
	mexPrintf("SYNTAX\n");
	mexPrintf("route=getroute(VER,ITRI,startpnt,endpnt)\n");
	mexPrintf("[route,length]=getroute(VER,ITRI,startpnt,endpnt)\n");
	mexPrintf("route=getroute(ADJ,startpnt,endpnt)\n");
	mexPrintf("[route,length]=getroute(ADJ,startpnt,endpnt)\n\n");
	mexPrintf("DESCRIPTION\n");		
	mexPrintf("route=getroute(VER,ITRI,startpnt,endpnt)  Returns the shortest route from startpnt to endpnt.\n");
	mexPrintf("                  startpnt and endpnt are vertex indices. The shortest route is based on the\n");
	mexPrintf("                  adjacency matrix weighed with the edge distances. If the graph is inconsitent \n");
	mexPrintf("                  and no route can be found route returns empty.\n");
	mexPrintf("[route,length]=getroute(VER,ITRI,startpnt,endpnt)  Returns the shortest route from startpnt to\n");
	mexPrintf("                  endpnt. Additionaly returns the length of this route. If the graph is inconsitent \n");
	mexPrintf("                  and no route can be found route and length return empty.\n");
	mexPrintf("route=getroute(ADJ,startpnt,endpnt)  Returns the shortest route from startpnt to endpnt.\n");
	mexPrintf("                  startpnt and endpnt are vertex indices of the ADJacency matrix. The shortest \n");
	mexPrintf("                  route is based on the distances given by the adjacency matrix. If the graph \n");
	mexPrintf("                  is inconsitent and no route can be found route returns empty.\n");
	mexPrintf("[route,length]=getroute(ADJ,startpnt,endpnt)  Returns the shortest route from startpnt to endpnt.\n");
	mexPrintf("                  Additionaly returns the length of this route If the graph is inconsitent \n");
	mexPrintf("                  and no route can be found route and length return empty.\n\n");
	mexPrintf("See also\n- graphdist\n");

}

	
/***************************************************************************/

int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
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
void createAdjacency2D(double *itri)
{
	int i,j,j1,j2;
	double r[3];
	for (t=0;t<ndhk;t++)
	{
		dhk[t][0]=(int)itri[t]-1;
		dhk[t][1]=(int)itri[t+ndhk]-1;		
		dhk[t][2]=(int)itri[t+2*ndhk]-1;	
		for (j=0;j<3;j++)
		{
			j1=icyc(j,3);j2=icyc(j+1,3);
			for (i=0;i<3;i++)
				r[i]=ver[dhk[t][j1]+i*npnt]-ver[dhk[t][j2]+i*npnt];
			D[dhk[t][j1]+dhk[t][j2]*npnt]=norm(r);
			D[dhk[t][j1]*npnt+dhk[t][j2]]=D[dhk[t][j1]+dhk[t][j2]*npnt];
		}
	}
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
			if (D[i*npnt+j]>0) 
				nrib++;

	rib   = mxCalloc(nrib,2*sizeof(int));   rib	  || outofmem();
	k=0;  maxlength = 0;
	/* create a list of all the ribs */
	for ( i=0;i<npnt;i++)
	 	for (j=i;j<npnt;j++)
    		if (D[i*npnt+j]>0)
			{
	 		 	rib[k][0]=i;
				rib[k][1]=j;
				maxlength = maxlength+D[i*npnt+j];
				k++;
			}
	fixed = mxCalloc(npnt,sizeof(char));    fixed || outofmem();
	for (i=0; i<npnt;i++)
	{	
		for (j=0;j<npnt;j++)
		{
			if (D[i+j*npnt]==0 )
				D[i+j*npnt]=maxlength;
		}
		fixed[i]=0;
	}
}
/***************************************************************************/

int CalculateRoute(void)
{
	int fixrib[MAXNODES],i,n,nfix,start,maxNodes=MAXNODES-1;
	
	init_getdistmat();
	frompnt = mxCalloc(npnt,sizeof(int));    frompnt || outofmem();
	start = startPnt;
	for (i=0; i<npnt;i++)
		fixed[i]=0;
	nfix=1;
	D[start+start*npnt]=0;
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
					D[rib[i][0]+start*npnt]+ D[rib[i][0]+rib[i][1]*npnt]<=length)
			{
				length = D[rib[i][0]+start*npnt]+ D[rib[i][0]+rib[i][1]*npnt];
				fixrib[n]=i;
				n++;
				}
			else if(!fixed[rib[i][0]] && fixed[rib[i][1]] &&
						D[rib[i][1]+start*npnt]+ D[rib[i][1]+rib[i][0]*npnt] <=length)
			{
				length = D[rib[i][1]+start*npnt]+ D[rib[i][1]+rib[i][0]*npnt];
				fixrib[n]=i;
				n++;
			}
		}
		if (n==0)
			return 0; //mexErrMsgTxt("Inconsistent graph (more than 1 part).");
		if (n>maxNodes)
			mexErrMsgTxt("Too many nodes (>=1000) with the same length.");
		for (i=0;i<n;i++)
		{
			
			if(fixed[rib[fixrib[i]][1]] &&
				 fabs(D[rib[fixrib[i]][1]+start*npnt]+ D[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < EPS)      
			{
				//mexPrintf("length %f  %d  %d \n",length,rib[fixrib[i]][1],rib[fixrib[i]][0]);
				fixed[rib[fixrib[i]][0]]=1;
				D[rib[fixrib[i]][0]+start*npnt]=length;
				D[start+rib[fixrib[i]][0]*npnt]=length;
				frompnt[rib[fixrib[i]][0]]=rib[fixrib[i]][1];
				if (rib[fixrib[i]][0]==endPnt)
					return 1;

			}
			else if (fabs(D[rib[fixrib[i]][0]+start*npnt]+ D[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < EPS)      
			{
//				mexPrintf("length %f  %d  %d \n",length,rib[fixrib[i]][0],rib[fixrib[i]][1]);
				fixed[rib[fixrib[i]][1]]=1;
				D[rib[fixrib[i]][1]+start*npnt]=length;
				D[start+rib[fixrib[i]][1]*npnt]=length;
				frompnt[rib[fixrib[i]][1]]=rib[fixrib[i]][0];
				if (rib[fixrib[i]][1]==endPnt)
					return 1;
			}
		}
		nfix=0;
		for (i=0;i<npnt;i++)
			if (fixed[i])
				nfix++;
	}
	return 1;
}


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  mrows,mcols,i,j,nFroms,found;
	double *itri,*dt,*len;

	/* Check input arguments */
	if(nrhs >4 || nrhs <3 ||nlhs >2) 
	{
		showHelp();
		if (nrhs==0)
			return;
		else
			mexErrMsgTxt("Check arguments");		
	}
	if (nrhs ==4)
	{
		// VER	
	 	mrows=mxGetM(prhs[0]);
		mcols=mxGetN(prhs[0]);  	 	
		if( mcols !=3)
			mexErrMsgTxt("size VER = [.., 3]");
		if(!(mxIsDouble(prhs[0]))) 
			mexErrMsgTxt("VER array must be of type double.");
		npnt = mrows;
		// ITRI
		mrows=mxGetM(prhs[1]);
	 	mcols=mxGetN(prhs[1]);  	 	
		if( mcols !=3)
			mexErrMsgTxt("size ITRI = [.., 3]");
		ndhk = mrows;
		itri = mxGetPr(prhs[1]);
		for (i=0;i<ndhk;i++)
			if (itri[i]-(int)itri[i]!=0 )
				mexErrMsgTxt("second parameter (ITRI) must contain index numbers only");
		
		// start / end
		if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )
			mexErrMsgTxt("start pnt must be an scalar");
		
		if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )
			mexErrMsgTxt("end pnt must be an scalar");

		ver  = mxGetPr(prhs[0]);

		startPnt=(int)mxGetScalar(prhs[2])-1;
		if (startPnt < 0  | startPnt >=npnt)
			mexErrMsgTxt("start pnt must be between 1 and number of vertices");
		endPnt=(int)mxGetScalar(prhs[3])-1;
		if (endPnt < 0  | endPnt >=npnt)
			mexErrMsgTxt("end pnt must be between 1 and number of vertices");

		/*Set output pointer to output matrix*/   
		D = mxCalloc(npnt*npnt,sizeof(double));   D	  || outofmem();
		dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();
		createAdjacency2D(itri);
		found=CalculateRoute();
		mxFree(dhk);
	}
	else
	{
		// ADJ
	 	mrows=mxGetM(prhs[0]);
		mcols=mxGetN(prhs[0]);  	 	
		if( mrows !=mcols)
			mexErrMsgTxt("size ADJ: nrows must be equal tro ncols");
		if(!(mxIsDouble(prhs[0]))) 
			mexErrMsgTxt("ADJ array must be of type double.");
		npnt = mrows;
		ndhk=0;
		// start / end
		if( mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1 )
			mexErrMsgTxt("start pnt must be an scalar");
		if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )
			mexErrMsgTxt("end pnt must be an scalar");
		startPnt=(int)mxGetScalar(prhs[1])-1;
		endPnt=(int)mxGetScalar(prhs[2])-1;
		if (startPnt < 0  || startPnt >=npnt)
			mexErrMsgTxt("start pnt must be between 1 and number of vertices");
		if (endPnt < 0  || endPnt >=npnt)
			mexErrMsgTxt("end pnt must be between 1 and number of vertices");

		dt= mxGetPr(prhs[0]);
		D = mxCalloc(npnt*npnt,sizeof(double));   D	  || outofmem();
		for (i=0; i<npnt;i++)
			for (j=0; j<npnt;j++)
				D[i+j*npnt]=dt[i+j*npnt];
		found=CalculateRoute();
	}
	mxFree(D);
	if (found)
	{
		nFroms=1;
		j=endPnt;
		while ( j!= startPnt )
		{
			j=frompnt[j];
			nFroms++;
		}
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(1,nFroms,mxREAL);
		D= mxGetPr(plhs[0]);		
		i=nFroms-1;
		j=endPnt;// hier gaat iets mis nu
		//ga van eind naar begin!! zet nog een return als het eindpunt gevonden is!!
		while ( j!= startPnt )
		{
			D[i--]=j+1;
			j=frompnt[j];
		}
		D[i]=startPnt+1;
		if(nlhs ==2) 
		{
			plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
			len= mxGetPr(plhs[1]);		
			len[0]=length;
		}
   }	
	else
	{
		plhs[0]=mxCreateDoubleMatrix(1,0,mxREAL);
		if(nlhs ==2) 
			plhs[1]=mxCreateDoubleMatrix(1,0,mxREAL);
	}
	mxFree(frompnt);
	mxFree(fixed);
	mxFree(rib);
}

/**************************** End of distgraph.c *****************************/


