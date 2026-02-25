/*************************** distgraph.c ***********************************/
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
#define MAXDHKVANRIB 20
#define MAXNODES 3000
/************************* Defintions and globals **************************/

int		nrib,npnt,ndhk,t;
double	maxlength;
double	length;
double (*ver)			= NULL;
double (*D)				= NULL;
double (*S)				= NULL;
char	 (*fixed)		= NULL;
int	 (*rib)[2]		= NULL;

/***************************************************************************/


void showHelp(void)
{
    mexPrintf("splitgraph \n determine for each node to which graph it belongs\n\n");
    mexPrintf("SYNTAX\n");
    mexPrintf("S=splitgraph (ADJ) \n");
    mexPrintf("S=splitgraph(VER,ITRI) \n");


    mexPrintf("DESCRIPTION\n");
    mexPrintf("S=splitgraph(VER,ITRI),  A triangulated surface can consist out of more parts. Each part is a separate graph\n");
    mexPrintf("                      for each of the nodes it is determined to which graph it belongs. \n");
    mexPrintf("                      S has the same length as the number of nodes. Each node get the number of the graph it belongs to.\n");
    mexPrintf("                      Input is the triangulated surface describetd by its vertices and vertex indexes (ITRI)\n");
    mexPrintf("S=splitgraph (ADJ), Same as previous, input is now the adjacency matrix \n");
}


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
    int	 (*dhk)[3]		=NULL;
    dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();

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
    mxFree(dhk);
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

int SplitGraph(void)
{
    int fixrib[MAXNODES],i,n,nfix,start,maxNodes=MAXNODES-1;
    int part=1;

    start = 1;
    init_getdistmat();
    for (i=0; i<npnt;i++)
        fixed[i]=0;
    nfix=1;
    fixed[start]=1;
    S[start] = part;
    D[start+start*npnt]=0;
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
        if (n>maxNodes)
            mexErrMsgTxt("Too many nodes (>=150) with the same length.");
        if (n==0)
        {
            for (i=0;i<npnt;i++)
                if (!fixed[i])
                {
                    start=i;
                    part=part+1;
                    S[start]=part;
                    fixed[start]=1;
                    D[start+start*npnt]=0;
                    break;
                }
            //mexErrMsgTxt("Inconsistent graph (more than 1 part).");
        }

        for (i=0;i<n;i++)
        {
            if(fixed[rib[fixrib[i]][1]] &&
                    fabs(D[rib[fixrib[i]][1]+start*npnt]+ D[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < EPS)
            {
                fixed[rib[fixrib[i]][0]]=1;
                D[rib[fixrib[i]][0]+start*npnt]=length;
                D[start+rib[fixrib[i]][0]*npnt]=length;
                S[rib[fixrib[i]][0]]=part;
            }
            else if (fabs(D[rib[fixrib[i]][0]+start*npnt]+ D[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < EPS)
            {
                fixed[rib[fixrib[i]][1]]=1;
                D[rib[fixrib[i]][1]+start*npnt]=length;
                D[start+rib[fixrib[i]][1]*npnt]=length;
                S[rib[fixrib[i]][1]]=part;
            }
        }
        nfix=0;
        for (i=0;i<npnt;i++)
            if (fixed[i])
                nfix++;
    }
    return part;
}


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Declare variables. */
    int  mrows,mcols,i,j,nFroms,part;
    double *itri,*dt,*len;

    /* Check input arguments */
    if (nrhs<1)
    {
        showHelp();
        return;
    }
    if(nrhs >2 )
    {
        showHelp();
        mexErrMsgTxt("Two input arguments required: 1) VER, 2) ITRI OR \n One input arguments required: 1) Adjacency matrix\n");
    }
    if(nlhs >1)
    {
        showHelp();
        mexErrMsgTxt("Too many output arguments. (one only)");
    }
    if (nrhs ==2)
    {
        // VER
        mrows=mxGetM(prhs[0]);
        mcols=mxGetN(prhs[0]);
        if( mcols !=3)
        {
            showHelp();
            mexErrMsgTxt("size VER = [.., 3]");
        }
        if(!(mxIsDouble(prhs[0])))
        {
            showHelp();
            mexErrMsgTxt("VER array must be of type double.");
        }
        npnt = mrows;
        // ITRI
        mrows=mxGetM(prhs[1]);
        mcols=mxGetN(prhs[1]);
        if( mcols !=3)
        {
            showHelp();
            mexErrMsgTxt("size ITRI = [.., 3]");
        }
        ndhk = mrows;
        ver  = mxGetPr(prhs[0]);
        itri = mxGetPr(prhs[1]);
        D = mxCalloc(npnt*npnt,sizeof(double));   D	  || outofmem();
        createAdjacency2D(itri);
        plhs[0]=mxCreateDoubleMatrix(1,npnt,mxREAL);
        S= mxGetPr(plhs[0]);
        part=SplitGraph();
        mxFree(D);
    }
    else
    {
        // ADJ
        mrows=mxGetM(prhs[0]);
        mcols=mxGetN(prhs[0]);
        if( mrows !=mcols)
        {
            showHelp();
            mexErrMsgTxt("size ADJ: nrows must be equal to ncols");
        }
        if(!(mxIsDouble(prhs[0])))
        {
            showHelp();
            mexErrMsgTxt("ADJ array must be of type double.");
        }
        npnt = mrows;
        dt= mxGetPr(prhs[0]);
        D = mxCalloc(npnt*npnt,sizeof(double));   D	  || outofmem();
        for (i=0;i<npnt;i++)
            for (j=0;j<npnt;j++)
                D[i+j*npnt]=dt[i+j*npnt];

        plhs[0]=mxCreateDoubleMatrix(1,npnt,mxREAL);
        S= mxGetPr(plhs[0]);
        part=SplitGraph();
    }
    mxFree(fixed);
    mxFree(rib);
    if (part<10)
        mexPrintf(" graph consists of %d parts\n",part);
}

/**************************** End of distgraph.c *****************************/


