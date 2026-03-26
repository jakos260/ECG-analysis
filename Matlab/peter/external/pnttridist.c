
  /*************************************************************************/
  /*									   */
  /*	pnttridist: program to compute distance between a point 	   */
  /*	and a triangluated surface.					   */
  /*									   */
  /*									   */
  /* 	20-10-00				Thom Oostendorp		   */
  /*						Medical Physics		   */
  /*						University of Nijmegen	   */
  /*									   */
  /*************************************************************************/



#include "mex.h"
#include <math.h>
//#include "C:\ECG_simulation\Thom\trilib.h"

#define DELTA 0.01
#define M_PI 3.14159265358979323846

double ddeter(double r1[3], double r2[3], double r3[3])
 {
  double det;
   
  det=r1[0]*(r2[1]*r3[2]-r2[2]*r3[1]);
  det=det-r1[1]*(r2[0]*r3[2]-r2[2]*r3[0]);
  det=det+r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]);
  return det;
 }
	 
double vecdist(double *r1, double *r2)
 {
  int	k;
  double r[3];
  
  for (k=0; k<3; k++)
      r[k]=r1[k]-r2[k];
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
 }
 
 	/*   plindist(p, l1, l2, cp, lab): distance beween point p and	*/
	/*	line segment l1-l2. In cp the point on l1-l2 closed to	*/
	/*	p is returned; in lab the fraction of the distance	*/
	/*	l1-cp in respect to l1-l2 is returned.			*/

double plindist(double p[3], double l1[3], double l2[3], double cp[3], double *lab1)
{
    int	k;
    double s[3];
    double sq=0, lab=0;
    
    for (k=0; k<3; k++)
    {
        s[k]=l2[k]-l1[k];
        sq += s[k]*s[k];
        lab -= (l1[k]-p[k])*s[k];
    }
    if (sq==0)
        lab=0;
    else
        lab /= sq;
    if (lab<0)
        lab=0;
    if (lab>1)
        lab=1;
    for (k=0; k<3; k++)
        cp[k]=l1[k]+lab*s[k];
    *lab1=lab;
    return vecdist(p,cp);
}
 
int lmoutr(double *p1, double *p2, double *p3, double *r, 
		double *lambda, double *mu, double *dif)
{
    double s1[3], s2[3], det;
    double n[3], bb[3], hulp[3], tmp; 
    int k, l;
    
    for (k=0; k<3; k++)
    {
        s1[k]=p2[k]-p1[k];
        s2[k]=p3[k]-p1[k];
        bb[k]=r[k]-p1[k];
    }
    for (k=0; k<3; k++)
    {
        for (l=0; l<3; l++)
        {
            hulp[l]=0;
        }
        hulp[k]=1;
        n[k]=ddeter(hulp,s1,s2);
    }
    tmp=n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
    if (tmp==0)
        return 0;
    tmp=sqrt(sqrt(tmp));
    n[0]=n[0]/tmp;
    n[1]=n[1]/tmp;
    n[2]=n[2]/tmp;
    det=ddeter(s1, s2, n);
    if (det==0)
        return 0;
    *lambda=ddeter(bb,s2,n)/det;
    *mu=ddeter(s1,bb,n)/det;
    *dif=ddeter(s1,s2,bb)/det;
    return 1;
}

int pntorder(double p1[3], double p2[3])
{
    int k;
    
    for (k=0; k<3; k++)
    {
        if (p1[k]>p2[k])
            return 1;
        else if (p1[k]<p2[k])
            return -1;
    }
    return 0;
}

double ptridist(double p[3], double v1[3], double v2[3], double v3[3],
		double cp[3], double *lab, double *mu)
{
    int	k;
    double d;
    
    lmoutr(v1, v2, v3, p, lab, mu, &d);
    if ( *lab>=0 && *mu>=0 && (*lab+*mu)<=1)
    {
        for (k=0; k<3; k++)
            cp[k]=(1-*lab-*mu)*v1[k]+*lab*v2[k]+*mu*v3[k];
        return vecdist(p, cp);
    }
    if (*lab<0)
    {
        *lab=0;
        return plindist(p, v1, v3, cp, mu);
    }
    if (*mu<0)
    {
        *mu=0;
        return plindist(p, v1, v2, cp, lab);
    }
    d=plindist(p, v2, v3, cp, mu);
    *lab=1-*mu;
    return d;
}



/* pntontotri projects r onto triangulated surface */

int pntontotri(double rin[3], double rpmin[3],int npnt, int ndhk, double (*pnt)[3], int (*dhk)[3],double *labMin,double *muMin)
{
	int	idhk, k, mindhk=0;
	double	dist, mindist=0, rp[3],lab,mu;
    
	mindist=ptridist(rin, pnt[dhk[0][0]], pnt[dhk[0][1]],pnt[dhk[0][2]], rpmin, &lab, &mu);
	for (idhk=1; idhk<ndhk; idhk++)
	{
		dist=ptridist(rin, pnt[dhk[idhk][0]], pnt[dhk[idhk][1]],pnt[dhk[idhk][2]], rp, &lab, &mu);
		if ( dist < mindist)
		{
			mindhk=idhk;
			mindist=dist;
			for (k=0; k<3; k++)
            {
				rpmin[k]=rp[k];
            }
            *labMin = lab;
            *muMin = mu; 
		}
	}
	return mindhk;
}


void print_help(void)
 {
  mexPrintf("");
  mexPrintf("  usage: pnttridist [dist,loc,onItri,LAMBDA,MU]=pnttridist(VER,ITRI,PNTS)\n");
  mexPrintf("");
  mexPrintf("  The minimal distances between the points in PNTS\n");
  mexPrintf("  and the surface described by VER and ITRI are computed \n\n");
  mexPrintf("  The output is put to dist and loc.\n");
}
int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}


/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{
  /* Declare variables. */  
	int  mrows1,mrows,mcols,mcols1,minIdhk;
	int 	i, j, t, npnt, npnt1, ndhk, (*dhk)[3];
	double (*pnt)[3], (*pnt1)[3];
    double pnt2[3];
    double lambdaP,muP;
	double *ver, *ver1, *DIST, *LOC, *itri, *VER,*LAMBDA,*MU;
  
	if (nrhs!=3)
	{
		print_help();
		mexErrMsgTxt("Inproper number of arguments ");
	}	
  	/* ITRI */ 
	mrows=mxGetM(prhs[1]);
 	mcols=mxGetN(prhs[1]);  	 	
	if( mcols !=3)						
		mexErrMsgTxt("size ITRI = [.., 3]");
	ndhk = mrows;

	/* VER	*/
 	mrows=mxGetM(prhs[0]);
	mcols=mxGetN(prhs[0]);  	 	
	if( mcols !=3)						
		mexErrMsgTxt("size VER = [.., 3]");
	if(!(mxIsDouble(prhs[0]))) 
		mexErrMsgTxt("VER array must be of type double.");
	npnt = mrows;

	mrows1=mxGetM(prhs[2]);
	mcols1=mxGetN(prhs[2]);  	 	
	if( mcols1 != 3)						
		mexErrMsgTxt("size PNTS = [.., 3]");
	if(!(mxIsDouble(prhs[2]))) 
		mexErrMsgTxt("VER array must be of type double.");
	npnt1 = mrows1;
	
	ver  = mxGetPr(prhs[0]);
	itri = mxGetPr(prhs[1]);
	ver1 = mxGetPr(prhs[2]);	
	for (i=0;i<ndhk;i++)
		if (itri[i]-(int)itri[i]!=0 )
			mexErrMsgTxt("first parameter (ITRI) must contain index numbers only");
	dhk=mxCalloc(ndhk, 3*sizeof(int));		dhk	|| outofmem();
	pnt=mxCalloc(npnt, 3*sizeof(double));	pnt	|| outofmem();
	pnt1=mxCalloc(npnt1, 3*sizeof(double));pnt1	|| outofmem();	
	for (t=0;t<ndhk;t++)
	{
		dhk[t][0]=(int)itri[t]-1;
		dhk[t][1]=(int)itri[t+ndhk]-1;		
		dhk[t][2]=(int)itri[t+2*ndhk]-1;	
	}

	for (t=0;t<npnt;t++)
	{
		pnt[t][0]=ver[t];
		pnt[t][1]=ver[t+npnt];		
		pnt[t][2]=ver[t+2*npnt];	
	}
	
	for (t=0;t<npnt1;t++)
	{
		pnt1[t][0]=ver1[t];
		pnt1[t][1]=ver1[t+npnt1];		
		pnt1[t][2]=ver1[t+2*npnt1];	
	}
	
	plhs[0]=mxCreateDoubleMatrix(npnt1,1,mxREAL);					
	DIST = mxGetPr(plhs[0]);
	if (nlhs>1)
	{
		plhs[1]=mxCreateDoubleMatrix(npnt1,3,mxREAL);		
		VER= mxGetPr(plhs[1]);
	}
	if (nlhs>2)
	{
		plhs[2]=mxCreateDoubleMatrix(npnt1,1,mxREAL);		
		LOC=mxGetPr(plhs[2]);
	}
    if (nlhs>3)
	{
		plhs[3]=mxCreateDoubleMatrix(npnt1,1,mxREAL);		
		LAMBDA=mxGetPr(plhs[3]);
	}
    if (nlhs>4)
	{
		plhs[4]=mxCreateDoubleMatrix(npnt1,1,mxREAL);		
		MU=mxGetPr(plhs[4]);
	}

	for (i=0; i<npnt1; i++)
	{
  		minIdhk=pntontotri(pnt1[i],pnt2,npnt, ndhk, pnt, dhk,&lambdaP,&muP);
		DIST[i]=vecdist(pnt1[i], pnt2);
		if (nlhs>1)
        {
			for (j=0;j<3;j++)
            {
				VER[i+j*npnt1]=pnt2[j];
            }
        }
		if (nlhs>2) 
        {
            LOC[i]=minIdhk+1;
        }
        if (nlhs>3) 
        {
//             mexPrintf("The %f \n",lambdaP);
            LAMBDA[i]=lambdaP;
        }
        if (nlhs>4) 
        {
//             mexPrintf("mu %f \n",muP);

            MU[i]=muP;            
        }
    }
	
	mxFree(pnt);
	mxFree(dhk);
	mxFree(pnt1);
}

