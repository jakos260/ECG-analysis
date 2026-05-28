/*************************** tri2graph.c ***********************************/
/*                                                                         */           
/* This program calculates the shortest paths between one vertex and all   */
/* the other vertices of a triangulated object.                            */
/* The program creates square matrix with the number of vertices.				*/
/* Input1: An adjaceny matrix (file) with the distances between the        */
/*         adjacent vertices                                               */
/* Input2: Vertices and the triangle indexes (VER, ITRI)                   */
/* Input3: mode determines the output
/*						1: Adjecency matrix with in wall connections
						2: Distace matrix based on mode 0 matrix 
						3: Distance matrix based on mode 1 matrix 
						4: Distace and paths matrix based on mode 0 matrix
						5: Distance and paths matrix based on mode 1 matrix 
/* Output:  A matrix with all shortest paths between the all vertices.     */
/*				A matrix with all the node use for everey calcylated focus		*/
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
#define MAXNODES 300
/************************* Defintions and globals **************************/

int		nrib,npnt,ndhk,t;
double	maxlength;
int	 (*dhk)[3]		=NULL;
double (*ver)			= NULL;
double (*D)				= NULL;
char	 (*fixed)		= NULL;
int	 (*rib)[2]		= NULL;
double (*vertex_use)	= NULL;

/***************************************************************************/

int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}

/***************************************************************************/

int compare( const void *elem1, const void *elem2 )
{
	return &elem1 <= &elem2;
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
double rhoek(double r1[3], double r2[3], double r3[3])
{
	double cp23_x,cp23_y,cp23_z,n1,n2,n3,ip12,ip23,ip13,nom,den;
	
	cp23_x = r2[1] * r3[2] - r2[2] * r3[1];
	cp23_y = r2[2] * r3[0] - r2[0] * r3[2];
	cp23_z = r2[0] * r3[1] - r2[1] * r3[0];
	nom = cp23_x*r1[0] + cp23_y*r1[1] + cp23_z*r1[2];
	n1 = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
	n2 = sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);
	n3 = sqrt(r3[0]*r3[0] + r3[1]*r3[1] + r3[2]*r3[2]);
	ip12 = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
	ip23 = r2[0]*r3[0] + r2[1]*r3[1] + r2[2]*r3[2];
	ip13 = r1[0]*r3[0] + r1[1]*r3[1] + r1[2]*r3[2];
	den = n1*n2*n3 + ip12*n3 + ip23*n1 + ip13*n2;
	if (nom==0)
		if (den<=0)
			return 0;
	return -2.*atan2 (nom,den);
}

/***************************************************************************/
double deter(double r1[3], double r2[3], double r3[3])
 {
  double det;
   
  det=r1[0]*(r2[1]*r3[2]-r2[2]*r3[1]);
  det=det-r1[1]*(r2[0]*r3[2]-r2[2]*r3[0]);
  det=det+r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]);
  return det;
  
}
/***************************************************************************/

int nodeview(int j1,int j2)
{
	/* view =1 if line connecting node j1 to node j2 of
	 the triangulated surface specified by (VER and ITRI) lies entirely in 
	 its interior; else view=0*/

	double a1[3],a2[3],a3[3],sav[3],lambda,mu;
	double r1[3],r2[3],r3[3],det,position[3],hoek;
	int i,n,na,ns,p; 
	double alphas[100], alph[100],alpha;

	
	for (i=0;i<3;i++)
		a3[i]=ver[j1+i*npnt]-ver[j2+i*npnt];
	ns=0;
	alpha=0;
	// identify intersections in between node1 and node2
	for (i=0;i<ndhk;i++)
	{
		if ( dhk[i][0]!=j1 || dhk[i][1]!=j1 || dhk[i][2]!=j1 || 
			  dhk[i][0]!=j2 || dhk[i][1]!=j2 || dhk[i][2]!=j2 )
		{	
			for (p=0;p<3;p++)
			{
				a1[p]=ver[dhk[i][1]+p*npnt]-ver[dhk[i][0]+p*npnt];
	 		   a2[p]=ver[dhk[i][2]+p*npnt]-ver[dhk[i][0]+p*npnt];
 			}
  	   	det=deter(a1, a2, a3);
	  	 	if (fabs(det) > 1.e-9)
	   	{
     	   	// no paralel detected
      	  	// solution via Cramer's rule
       	 	for (p=0;p<3;p++)
	      		sav[p]=ver[j1+p*npnt]-ver[dhk[i][0]+p*npnt];
		      lambda=deter(sav,a2,a3)/det;
 	       	if (lambda >= 0 && lambda <= 1)
  	      	{
   	        	mu=deter(a1, sav, a3)/det;
    	       	if (mu >= 0 && mu <= 1 && (lambda + mu) <= 1)
     	      	{
      	      	alpha=deter(a1,a2,sav)/det; 
       	       	if (alpha > 0 && alpha < 1)
        	      	{
         	      	ns=ns+1; 
          	      	alphas[ns]=alpha;            
           	 	   }
					}	
				}
			}
		}
	}
	
	ns=ns+1;
	alphas[ns]=1;

	// sort for increasing values of alpha
	qsort( alphas, ns, sizeof(int), compare );
	na=0;
	alph[na++]=alphas[0];
	for (n=1;n<ns;n++)
		if (alphas[n]!=alphas[n-1])
			alph[na++]=alphas[n]; 
	alph[na++]=1;
	// test if positions halfway the intersections are interior or exterior
	// if any of these is exterior, nodes j1 and j2 cannot be connected
	
	for (n=1; n < na;n++)
	{
		for(p=0;p<3;p++)
			position[p]=ver[j1+p*npnt]+((alph[n]+alph[n-1])/2)*(ver[j2+p*npnt]-ver[j1+p*npnt]);//     a3[p];
		hoek=0;	
		for (i=0;i<ndhk;i++)
		{
			for(p=0;p<3;p++)
			{
				r1[p]=ver[dhk[i][0]+npnt*p]-position[p];
				r2[p]=ver[dhk[i][1]+npnt*p]-position[p];
				r3[p]=ver[dhk[i][2]+npnt*p]-position[p];			
			}	
			hoek+=rhoek(r1,r2,r3);
		}
		if (fabs(hoek) < 0.1)
			return 0; 
	}
	return 1;
}

/***************************************************************************/
void createAdjacency()
{
	int j,j1,j2;
	for (t=0;t<ndhk;t++)
		for (j=0;j<3;j++)
		{
			j1=icyc(j,3);j2=icyc(j+1,3);
			D[dhk[t][j1]+dhk[t][j2]*npnt]=1;
			D[dhk[t][j1]*npnt+dhk[t][j2]]=1;
		}
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
void calcAdjacency3D(void)
{
	int j1,j2,i,view;
	double r[3];
	for (j1=0;j1<npnt;j1++)
	{
		mexPrintf("node %d \n",j1);
		for (j2=j1+1; j2<npnt;j2++)
		{
			view=nodeview(j1,j2);
			if (view==1)
			{
				for (i=0;i<3;i++)
					r[i]=ver[j1+i*npnt]-ver[j2+i*npnt];
				D[j1+npnt*j2]=norm(r);
				D[j2+npnt*j1]=D[j1+npnt*j2];
			}
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

void CalculateDists(void)
{
	double length;
	int fixrib[MAXNODES],i,n,nfix,start, maxNodes=MAXNODES-1;
	init_getdistmat();
	for (start=0; start<npnt;start++)
	{
		for (i=0; i<npnt;i++)
			fixed[i]=0;
		nfix=1;
		D[start+start*npnt]=0;
		fixed[start]=1;
		/* calculate all paths between startpnt and vertices*/
		while (nfix<npnt)
		{
			length = maxlength-0.01;
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
				mexErrMsgTxt("Inconsistent graph (more than 1 part).");
			if (n>maxNodes)
				mexErrMsgTxt("Too many nodes (>=150) with the same length.");
			for (i=0;i<n;i++)
			{
				if( fixed[rib[fixrib[i]][1]] &&
					 fabs(D[rib[fixrib[i]][1]+start*npnt]+ D[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < EPS)      
				{
					fixed[rib[fixrib[i]][0]]=1;
					D[rib[fixrib[i]][0]+start*npnt]=length;
					D[start+rib[fixrib[i]][0]*npnt]=length;
					
					
				}
				else if (fabs(D[rib[fixrib[i]][0]+start*npnt]+ D[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < EPS)      
				{
					fixed[rib[fixrib[i]][1]]=1;
					D[rib[fixrib[i]][1]+start*npnt]=length;
					D[start+rib[fixrib[i]][1]*npnt]=length;
				}
			}
			nfix=0;
			for (i=0;i<npnt;i++)
				if (fixed[i]) 
					nfix++;
		}
	}
}
/***************************************************************************/

void CalculateDistsAndPaths(void)
{
	double length;
	int fixrib[MAXNODES],i,j,n,nfix,start,maxNodes=MAXNODES-1;
	int  (*frompnt)=NULL;
	
	init_getdistmat();
	frompnt = mxCalloc(npnt,sizeof(int));    frompnt || outofmem();
	for (start=0; start<npnt;start++)
	{
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
				mexErrMsgTxt("Inconsistent graph (more than 1 part).");
			if (n>maxNodes)
				mexErrMsgTxt("Too many nodes (>=150) with the same length.");
			for (i=0;i<n;i++)
			{
				if(fixed[rib[fixrib[i]][1]] &&
					 fabs(D[rib[fixrib[i]][1]+start*npnt]+ D[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < EPS)      
				{
					fixed[rib[fixrib[i]][0]]=1;
					D[rib[fixrib[i]][0]+start*npnt]=length;
					D[start+rib[fixrib[i]][0]*npnt]=length;
					frompnt[rib[fixrib[i]][0]]=rib[fixrib[i]][1];
				}
				else if (fabs(D[rib[fixrib[i]][0]+start*npnt]+ D[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < EPS)      
				{
					fixed[rib[fixrib[i]][1]]=1;
					D[rib[fixrib[i]][1]+start*npnt]=length;
					D[start+rib[fixrib[i]][1]*npnt]=length;
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
			while ( j!= start )
			{
				vertex_use[j+start*npnt]++;
				j=frompnt[j];
			}
		}
		vertex_use[start+start*npnt]=npnt;
	}
	mxFree(frompnt);
}


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  mrows,mcols,mode,i,j;
	int j1,j2;
	double *itri,*dt;
	
	/* Check input arguments */
	if(nrhs >3 || nrhs <1 ) 
	{
		mexPrintf("Tri2graph calculates the shortest path in a fully connected graph\n");
		mexPrintf("The following 6 options are available: \n");
		mexPrintf(" - mode 0         : Adjacency matrix weighted by the eddge distance \n");
		mexPrintf(" - mode 1         : Adjacency type of matrix: the distances used being those between individual nodes, provided \n                    that they can ''see'' one another in the interior of  the triangulated surface \n");
		mexPrintf(" - mode 2         : Distance matrix based on the mode 0 type adjacency matrix \n");
		mexPrintf(" - mode 3         : Distance matrix based on the mode 1 type adjacency matrix \n");
		mexPrintf(" - mode 4         : Distance matrix as in mode 2; extra: shortest paths between individual nodes \n");
		mexPrintf(" - mode 5         : Distance matrix as in mode 3; extra: shortest paths between individual nodes \n");
		mexPrintf(" - mode undefined : Adjacency matrix. One Input required only ( adjacency=tri2graph(ITRI) )\n\n");
  		mexErrMsgTxt("Three input arguments required: 1) VER, 2) ITRI 3) mode   OR \n    Two input arguments required: 1) ADJ, 2) modeOR \n    One input argumet required 1) ITRI \n\n\n");
	}
  	if(nlhs >2 || (nrhs ==1 && nlhs>1) ) 
	{
		mexPrintf("Tri2graph calculates the shortest path in a fully connected graph\n");
		mexPrintf("The following 6 options are available: \n");
		mexPrintf(" - mode 0         : Adjacency matrix weighted by the eddge distance \n");
		mexPrintf(" - mode 1         : Adjacency type of matrix: the distances used being those between individual nodes, provided \n                    that they can ''see'' one another in the interior of  the triangulated surface \n");
		mexPrintf(" - mode 2         : Distance matrix based on the mode 0 type adjacency matrix \n");
		mexPrintf(" - mode 3         : Distance matrix based on the mode 1 type adjacency matrix \n");
		mexPrintf(" - mode 4         : Distance matrix as in mode 2; extra: shortest paths between individual nodes \n");
		mexPrintf(" - mode 5         : Distance matrix as in mode 3; extra: shortest paths between individual nodes \n");
		mexPrintf(" - mode undefined : Adjacency matrix. One Input required only ( adjacency=tri2graph(ITRI) )\n\n");		
		mexErrMsgTxt("Too many output arguments.\n\n\n");
	}
	if (nrhs ==3)
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
		// mode
		if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )
			mexErrMsgTxt("mode must be an scalar");
		
		ver  = mxGetPr(prhs[0]);
		itri = mxGetPr(prhs[1]);
		mode=(int)mxGetScalar(prhs[2]);
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);
		D = mxGetPr(plhs[0]);
		dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();

		createAdjacency2D(itri);

		if (mode==1 || mode==3 || mode==5 )
			calcAdjacency3D();	
	}
	else if (nrhs ==2)
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
		// mode
		if( mxGetM(prhs[1])!=1 && mxGetN(prhs[1])!=1 )
			mexErrMsgTxt("mode must be an scalar");
		dt  = mxGetPr(prhs[0]);
		ver  = NULL;
		itri = NULL;
		mode=(int)mxGetScalar(prhs[2]);
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);
		D = mxGetPr(plhs[0]);
		for (i=0; i<npnt;i++)
			for (j=0; j<npnt;j++)
				D[i+j*npnt]=dt[i+j*npnt];
	}
	else //nrhs ==1
	{
		mode=-1;
		// ITRI
		mrows=mxGetM(prhs[0]);
	 	mcols=mxGetN(prhs[0]);  	 	
		if( mcols !=3)
			mexErrMsgTxt("size ITRI = [.., 3]");
		ndhk = mrows;
		dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();
		itri = mxGetPr(prhs[0]);
		npnt=0;
		for (t=0;t<ndhk;t++)
		{
			dhk[t][0]=(int)itri[t]-1;
			dhk[t][1]=(int)itri[t+ndhk]-1;		
			dhk[t][2]=(int)itri[t+2*ndhk]-1;	
			if (dhk[t][0]>npnt) npnt=dhk[t][0];
			if (dhk[t][1]>npnt) npnt=dhk[t][1];
			if (dhk[t][2]>npnt) npnt=dhk[t][2];
		}
		npnt++;
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);
		D = mxGetPr(plhs[0]);

		createAdjacency();
		rib=NULL;
		fixed=NULL;
	}
	if ((mode==4 || mode==5) && nlhs==2 )
	{
		mexPrintf("Calculate distances and paths\n");
		plhs[1] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);
		vertex_use = mxGetPr(plhs[1]);
		CalculateDistsAndPaths();
	}
	else if (mode==2 || mode==3 || ((mode==4 || mode==5) && nlhs==1 ))
	{
		mexPrintf("Calculate distances\n");
		CalculateDists();
	}	
	mxFree(fixed);
	mxFree(rib);
	mxFree(dhk);
}

/**************************** End of tri2graph.c *****************************/

