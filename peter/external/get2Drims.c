/****************************** get2Drimbs.c ********************************/
/*                                                                         */           
/* This program calculates the shortest paths between one vertex and all   */
/* the other vertices of a triangulated object.                            */
/* The program creates matrix of 1 by number of vertices.									 */
/* Input:  An adjaceny matrix (file) with the distances between the        */
/*         adjacent vertices                                               */
/* Output: A matrix with all shortest paths between the all vertices.      */
/*                                                                         */
/* 09-09-03   Peter van Dam                                                */
/*            Medical Physics		                                           */
/*			      University of Nijmegen	                                     */  
/*                                                                         */
/***************************************************************************/

#include <mex.h>

/************************* Defintions and globals **************************/

int  nrib,npnt;
int  (*rib)[2] = NULL;
int  (*rimlst)=NULL;

/***************************************************************************/

int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}

/***************************************************************************/
void CreateRiblst(const mxArray *array_ptr)
{
  int i,j,k,found,tot;
  double *pr;
	int *ir,*jc,start_row,stop_row,row,col,n,tnrib;
	int  (*trib)[2] = NULL;
	
 	nrib=0;
 	if (mxIsSparse(array_ptr))
 	{
 		pr=mxGetPr(array_ptr);
 		ir=mxGetIr(array_ptr);
 		jc=mxGetJc(array_ptr);
 		n =mxGetN(array_ptr);
		/* create a list of all the ribs */
		tnrib = mxGetNzmax(array_ptr);
		nrib = tnrib/2;
  	rib  = mxCalloc(nrib,2*sizeof(int));   rib 	|| outofmem();
		k=0;  tot=0;
		for (col=0; col<n; col++)  
 		{ 
   		start_row = jc[col];   
   		stop_row  = jc[col+1]; 
     	for (row = start_row; row < stop_row; row++)  
			{
				found =0;
				for (j=0;j<k;j++)
					if ( (rib[j][0]==col && rib[j][1]==ir[row]) || 
							 (rib[j][0]==ir[row] && rib[j][1]==col) )
						found=1;
				if (!found)
				{
   		   	rib[k][0]=col;
					rib[k++][1]=ir[row];
		    	if (k>nrib)
		     	 	mexErrMsgTxt("Input matrix must be symetric (more).");
				}
				tot++;
			}
		}
		if (k!=nrib)
			mexErrMsgTxt("Input matrix must be symetric(less).");
 	}
 	else /* full matrix*/
 	{
 	  if(!(mxIsDouble(array_ptr))) 
  		mexErrMsgTxt("Input array must be of type double.");
 		pr = mxGetPr(array_ptr);
		for ( i=0;i<npnt;i++)
 		 	for (j=i;j<npnt;j++)
				if (pr[i*npnt+j]>0) 
					nrib++;
					
		rib  = mxCalloc(nrib,2*sizeof(int));   rib 	|| outofmem();
		k=0;
		/* create a list of all the ribs */
		for ( i=0;i<npnt;i++)
 		 	for (j=i;j<npnt;j++)
     		if (pr[i*npnt+j]>0)
				{
		 		 	rib[k][0]=i;
					rib[k][1]=j;
					k++;
				}
	}		 	
}

/***************************************************************************/
int findrib(int node1,int node2)
{
  int p;
  for(p=0;p<nrib;p++) 
  {
    if ( (rib[p][0]==node1 && rib[p][1]==node2) ||
    		 (rib[p][1]==node1 && rib[p][0]==node2) )
      return p;
  }
  return -1;
}

/***************************************************************************/
void CreateRimlst(const mxArray *array_ptr)
{
  int i,k,rows;
  double *pr;
	
 	pr=mxGetPr(array_ptr);
 	rows=mxGetM(array_ptr);
  rimlst=mxCalloc(nrib,sizeof(int));   rimlst 	|| outofmem();
	for (i=0;i<rows;i++)
	{
		k=findrib(pr[i]-1,pr[i+rows]-1);
		if (k >=0)
			rimlst[k]++;
		k=findrib(pr[i]-1,pr[i+rows+rows]-1);
		if (k >=0)
			rimlst[k]++;
		k=findrib(pr[i+rows]-1,pr[i+rows+rows]-1);
		if (k >=0)
			rimlst[k]++;
	}
	for (i=0; i<nrib;i++)
		if (rimlst[i] > 1)
			rimlst[i]=0;
}

/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{ 
  /* Declare variables. */  
  int  number_of_dims,mrows,mcols,i;
  double *outputptr;
    
  /* Check for proper number of input and output arguments. */    
  if(nrhs <1 ) 
  		mexErrMsgTxt("Two input arguments required ( 1) adjacency matrix, 2) Triagles matrix).");
  if(nlhs >2) 
   	mexErrMsgTxt("Too many output arguments.");
  /* Check data type of input argument. */

 	number_of_dims = mxGetNumberOfDimensions(prhs[0]);
 	if( number_of_dims != 2)
  	mexErrMsgTxt("First input matrix must be a 2 dimension matrix");

  mrows=mxGetM(prhs[0]);
  mcols=mxGetN(prhs[0]);  	 	
  if( mrows!=mcols)
  	mexErrMsgTxt("Input array must be a square matrix");

  npnt = mrows;
    		
  /* initialise */ 
  CreateRiblst(prhs[0]);  
	/*Set output pointer to output matrix*/   
  plhs[0]=mxCreateDoubleMatrix(nrib,2,mxREAL); 
  outputptr = mxGetPr(plhs[0]);
  for (i=0;i<nrib;i++)
  {
  	outputptr[i]=rib[i][0]+1;
  	outputptr[i+nrib]=rib[i][1]+1;  	
  }

	/* Calculate the rims */	
	if (nlhs ==2)
	{
		CreateRimlst(prhs[1]); 
	  plhs[1]=mxCreateDoubleMatrix(nrib,1,mxREAL);
	  outputptr = mxGetPr(plhs[1]);	  
	  for (i=0;i<nrib;i++)
  	  outputptr[i]=rimlst[i];
	}
	mxFree(rib);
	mxFree(rimlst);
}

/**************************** End of get2DRimbs.c ***************************/
