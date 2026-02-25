/*************************** distgraph.c ***********************************/
/*                                                                         */           
/* linetrisect(p1, p2, r1, r2, r3, alfa, lab, mu, parallel):					*/
/*	Detect wether the line segment between p1 and p2 and							*/
/*	the triangle r1, r2, r3 intersect. If so, the fraction						*/
/*	from p1 to p2 of the section point is returned as alfa,						*/
/*	and the fractions from r1 to r2 and r3 are returned							*/
/*	as lab and mu. In any case, parallel is set.										*/
/*	If either the line is degenerate, the function returns						*/
/*	-1; if the triangle is, it returns -2.												*/
/* Based on Software from Thom Oostendorp												*/
/* 10-05-05   Peter van Dam                                                */
/*                                                                         */
/***************************************************************************/

#include <mex.h>
#include <math.h>

typedef double VER[3];
typedef int		ITRI[3];


int outofmem(void)
{
  mexErrMsgTxt("Out of memory.");
  return 1;
}


double deps(void)
 {
  volatile double e, d;
  
  for (e=1.0, d=2.0; d>1.0; e /= 2)
    d=e+1;
  return 2*e;
 }
 
 double ddeter(double r1[3], double r2[3], double r3[3])
 {
  double det;
   
  det=r1[0]*(r2[1]*r3[2]-r2[2]*r3[1]);
  det=det-r1[1]*(r2[0]*r3[2]-r2[2]*r3[0]);
  det=det+r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]);
  return det;
 }
	
 void dveccross(double *r1, double *r2, double *r3)
 {
  r3[0]=r1[1]*r2[2]-r1[2]*r2[1];
  r3[1]=r1[2]*r2[0]-r1[0]*r2[2];
  r3[2]=r1[0]*r2[1]-r1[1]*r2[0];
 }
 
	/*	vecnorm(r):		vector length			*/
	
double dvecnorm(double *r)
 {
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
 }

int dsamevec(double r1[3], double r2[3])
 {
  int k;

  for (k=0; k<3; k++)
    if (r1[k]!=r2[k])
      return 0;
  return 1;
 }

 int dlinesect(double *p1, double *r1, double *p2, double *r2, double *lab,
		double *mu)
 {
  int	c1=0, c2=1;
  double det, z[3], a1, a2;

  dveccross(r1, r2, z);
  a1=z[2]*z[2];
  a2=z[1]*z[1];
  if (a2>a1)
   {
    a1=a2;
    c2=2;
   }
  a2=z[0]*z[0];
  if (a2>a1)
   {
    c1=1;
    c2=2;
   }

  det=r1[c2]*r2[c1]-r1[c1]*r2[c2];
  if (det==0)
    return 0;
  *lab=(r2[c2]*(p1[c1]-p2[c1])-r2[c1]*(p1[c2]-p2[c2]))/det;
  *mu =(r1[c2]*(p1[c1]-p2[c1])-r1[c1]*(p1[c2]-p2[c2]))/det;
  return 1;
 }

 
int dlmoutr(double *p1, double *p2, double *p3, double *r, 
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

int dlinetrisect(double p1i[3], double p2i[3], 
		 double r1i[3], double r2i[3], double r3i[3], 
		 double *alfa, double *lab, double *mu, int *parallel,
		 double p[3])
 {
  int	i, k, fs, pntswap;
  double p1[3], p2[3], diff;
  double r123[3][3];
  double *r1=r123[0], *r2=r123[1], *r3=r123[2];
  double det, dettmp, r1mp1[3], r[4][3], d[4][3], a, l, m;
  static int first=1;
  static double meps;

  if (first)
   {
    first=0;
    meps=100*deps();
   }
  pntswap=(pntorder(p1i, p2i)==1);
  for (k=0; k<3; k++)
   {
    if (pntswap)
     {
      p1[k]=p2i[k];
      p2[k]=p1i[k];
     }
    else
     {
      p1[k]=p1i[k];
      p2[k]=p2i[k];
     }
    r1[k]=r1i[k];
    r2[k]=r2i[k];
    r3[k]=r3i[k];
   }
  qsort(r123, 3, 3*sizeof(double), (int (*)(const void *, const void *))pntorder);
  for (k=0; k<3; k++)
   {
    r1mp1[k]=r1[k]-p1[k];
    d[0][k]=p2[k]-p1[k];
    d[1][k]=r1[k]-r2[k];
    d[2][k]=r1[k]-r3[k];
   }
  dveccross(d[1], d[2], d[3]);			
  if ((l=dvecnorm(d[3]))==0) return -1; 		/* triangle area 0 */
  l=dvecnorm(d[0])/l;
  det=ddeter(d[0], d[1], d[2]);
  for (k=0; k<3; k++)
    d[3][k]*=l;
  dettmp=ddeter(d[3], d[1], d[2]);
  if (dettmp==0)
   {
	mexPrintf("%8.4f %8.4f %8.4f\n", p1[0], p1[1], p1[2]);
	mexPrintf("%8.4f %8.4f %8.4f\n", p2[0], p2[1], p2[2]);
	mexPrintf("%8.4f %8.4f %8.4f\n", r1[0], r1[1], r1[2]);	
	mexPrintf("%8.4f %8.4f %8.4f\n", r2[0], r2[1], r2[2]);
	mexPrintf("%8.4f %8.4f %8.4f\n", r3[0], r3[1], r3[2]);
    mexErrMsgTxt("Special case in linetrisect; not yet implemented\n");

   }
  if (fabs(det/dettmp)<meps)
    det=0;
  if (det==0)
   {
    *parallel=1;
    if (dsamevec(p1, p2)) return -1;		/* line length 0 */
    dlmoutr(r1, r2, r3, p1, &l, &m, &a);
    if (a!=0)
      return 0;					/* l and t not in same plane */
    if (l>=0 && m>=0 && l+m <=1)
      return 1;					/* p1 in triangle */
    dlmoutr(r1, r2, r3, p2, &l, &m, &a);
    if (l>=0 && m>=0 && l+m <=1)
      return 1;					/* p2 in triangle */
    for (k=0; k<3; k++)
     {
      r[1][k]=r2[k];
      r[2][k]=r3[k];
      r[3][k]=r3[k];
      d[3][k]=r2[k]-r3[3];
     }
    for (i=1; i<3; i++)
     {
      if (!dlinesect(p1, d[0], r[i], d[i], &l, &m))	/* l and edge par. */
       {
	if (d[i][0]!=0)
	  fs=0;
	else if (d[i][1]!=0)
	  fs=1;
	else
	  fs=2;
        l=(p1[fs]-r[i][fs])/d[i][fs];
	for (k=fs+1; k<3; k++)
	  if (l!=(p1[k]-r[i][k])/d[i][k]);
	    goto nosect;
        m=(p2[fs]-r[i][fs])/d[i][fs];
	if ( (l<0 && m>1) || (l>1 && m<0) )	/* p1 and p2 at both sides */
	  return 1;
       }
      else
        if (l>0 && l<1 && m>0 && m<1)
	  return 1;				/* line crosses edge */
nosect: ;
     }
    return 0;					/* same plane, but no cross */
   }
  else
    *parallel=0;
  a=ddeter(r1mp1, d[1], d[2])/det;
  l=ddeter(d[0], r1mp1, d[2])/det;
  m=ddeter(d[0], d[1], r1mp1)/det;
  if (p)
    for (k=0; k<3; k++)
      p[k]=(1-a)*p1[k]+a*p2[k];
  //printf("%23.20f %23.20f %23.20f %23.30f \n\n", a, p[0], p[1], p[2]);
  if (pntswap)
    *alfa=1-a;
  else
    *alfa=a;
  //*lab=l;
  //*mu=m;
  dlmoutr(r1i, r2i, r3i, p, lab, mu, &diff);  // don't use l,m as vertices may be interchanged
  if ( a<0 || a>1 || l<0 || m<0 || l+m>1 )
    return 0;
  return 1;
 }

int linetrisect(double p1[3], double p2[3], 
		double r1[3], double r2[3], double r3[3], 
		double *alfa, double *lab, double *mu, int *parallel,
		double p[3])
 {
  int k, ret;
  double dp1[3], dr1[3], dp2[3], dr2[3], dr3[3], dl, dm, da, dp[3];

  for (k=0; k<3; k++)
   {
    dp1[k]=p1[k];
    dr1[k]=r1[k];
    dp2[k]=p2[k];
    dr2[k]=r2[k];
    dr3[k]=r3[k];
   }
  ret=dlinetrisect(dp1, dp2, dr1, dr2, dr3, &da, &dl, &dm, parallel, dp);
  *alfa=(double) da;
  *lab=(double) dl;
  *mu=(double) dm;
  if (p)
    for (k=0; k<3; k++)
      p[k]=(double) dp[k];
  return ret;
 }


/***************************************************************************/

/* mexFunction is the gateway routine for the MEX-file. */ 
void mexFunction( int nlhs, mxArray *plhs[],
             		int nrhs, const mxArray *prhs[] )
{

  /* Declare variables. */  
	int  t,mrows,mcols,i,j,n,jdhk,nerrors,npnt,ndhk,*orient,*dhki,sfound,efound;
	double *itri,*ver,*mu,*alpha,*lamb,meps,startPnt[3],endPnt[3],*td,*TRINTER;
	VER *pnt;
	ITRI *dhk;
		
	/* Check input arguments */
	if(nrhs !=4 ) 
	{
		mexPrintf("Inputs: \n - Vertices (nx3)\n - Triangle indices (nx3)\n - Two vertices describing the line (l1 and l2) \n\n");
		mexPrintf(" Outputs: TRIINTER\n");
		mexPrintf(" - TRIINTER(:,1): (index)  triangle indexes of intersection \n");
		mexPrintf(" - TRIINTER(:,2): (view)   parallel=0, view from outside=1, else -1 \n");
		mexPrintf(" - TRIINTER(:,3): (lambda) fraction along edge 1\n");
		mexPrintf(" - TRIINTER(:,4): (mu)     fraction along edge 2\n");
		mexPrintf(" - TRIINTER(:,5): (alpha)  fraction along the line from l1 to l2\n\n");
  		mexErrMsgTxt("Four input arguments required: 1) VER, 2) ITRI 3) first point  4) second point\n");
	}
  	if(nlhs >1) 
	{
		mexPrintf("Inputs: \n - Vertices (nx3)\n - Triangle indices (nx3)\n - Two vertices describing the line (l1 and l2) \n\n");
		mexPrintf(" Outputs: TRIINTER\n");
		mexPrintf(" - TRIINTER(:,1): (index)  triangle indexes of intersection \n");
		mexPrintf(" - TRIINTER(:,2): (view)   parallel=0, view from outside=1, else -1 \n");
		mexPrintf(" - TRIINTER(:,3): (lambda) fraction along edge 1\n");
		mexPrintf(" - TRIINTER(:,4): (mu)     fraction along edge 2\n");
		mexPrintf(" - TRIINTER(:,5): (alpha)  fraction along the line from l1 to l2\n\n");		
   	mexErrMsgTxt("Too many output arguments.");
	}
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
	// start / end
	if( mxGetM(prhs[2])==1 && mxGetN(prhs[2])==3  )
		td=mxGetPr(prhs[2]);
	else if ( mxGetM(prhs[2])==3 && mxGetN(prhs[2])==1  )
		td=mxGetPr(prhs[2]);
	else
		mexErrMsgTxt("start pnt must be an 3d vector");
		
	for(t=0;t<3;t++)
		startPnt[t]=td[t];
	
	if( mxGetM(prhs[3])==1 && mxGetN(prhs[3])==3 )
		td=mxGetPr(prhs[3]);
	else if( mxGetM(prhs[3])==3 && mxGetN(prhs[3])==1 )
		td=mxGetPr(prhs[3]);
	else
		mexErrMsgTxt("end pnt must be an 3d vector");
	for(t=0;t<3;t++)
		endPnt[t]=td[t];
	
	ver  = mxGetPr(prhs[0]);
	itri = mxGetPr(prhs[1]);
	dhk=mxCalloc(ndhk, 3*sizeof(int));		dhk	|| outofmem();
	pnt=mxCalloc(npnt, 3*sizeof(double));	pnt	|| outofmem();
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

	alpha=	mxCalloc(ndhk,sizeof(double));	alpha || outofmem();
	mu=		mxCalloc(ndhk,sizeof(double));	mu		|| outofmem();
	lamb=		mxCalloc(ndhk,sizeof(double));	lamb	|| outofmem();
	dhki=		mxCalloc(ndhk,sizeof(int));		dhki	|| outofmem();
	orient=	mxCalloc(ndhk,sizeof(int));		orient|| outofmem();
	for (jdhk=0; jdhk<ndhk; jdhk++)
	{
		alpha[jdhk]=0;
		mu[jdhk]=0;
		lamb[jdhk]=0;
		dhki[jdhk]=0;
		orient[jdhk]=0;
	}
	nerrors=0;
	for (jdhk=0; jdhk<ndhk; jdhk++)
	{			
		n=0;
		sfound=0;
		efound=0;
		for (i=0; i<3; i++)
			if((startPnt[0]==pnt[dhk[jdhk][i]][0] && 
				 startPnt[1]==pnt[dhk[jdhk][i]][1] &&
				 startPnt[2]==pnt[dhk[jdhk][i]][2]) )
			{
				sfound=1;
				break;
			}
		for (i=0; i<3; i++)		
			if((endPnt[0]  ==pnt[dhk[jdhk][i]][0] && 
				 endPnt[1]  ==pnt[dhk[jdhk][i]][1] &&
				 endPnt[2]  ==pnt[dhk[jdhk][i]][2]))
			{
				efound=1;
				break;
			}
		if (sfound || efound) 	continue;	/* end(s) of edge coincides with vertex(es) */

		i=linetrisect(startPnt, endPnt,
						  pnt[dhk[jdhk][0]], pnt[dhk[jdhk][1]], pnt[dhk[jdhk][2]], 
							&alpha[jdhk], &lamb[jdhk], &mu[jdhk], &orient[jdhk],NULL);
		if (i==-1)
			continue;	/* zero length edge; will be detected as zero area tr.*/
		if (i==-2)
			continue;	//Zero area triangle %4d\n", jdhk
		else if (i==1)
		{        
						// "The edge between vertices %4d and %4d intersects triangle jdhk
			dhki[jdhk]=jdhk;
			nerrors++;
		}
	}
	/*Set output pointer to output matrix*/   
	plhs[0]=mxCreateDoubleMatrix(nerrors,5,mxREAL);
	TRINTER= mxGetPr(plhs[0]);
	i=0;
	meps=1000*deps();
	for (jdhk=0; jdhk<ndhk; jdhk++)
	{	
		if (dhki[jdhk] !=0)
		{
			TRINTER[i+0*nerrors]=dhki[jdhk]+1;
			switch (orient[jdhk])
			{
				case 0 :{TRINTER[i+1*nerrors]=0;  break;}
				case -1:{TRINTER[i+1*nerrors]=1;  break;}
				default:{TRINTER[i+1*nerrors]=-1; break;}
			}
			TRINTER[i+1*nerrors]=orient[jdhk];
			TRINTER[i+2*nerrors]=lamb[jdhk]<meps ? 0 : fabs(lamb[jdhk]-1) <meps ? 1:lamb[jdhk];
			TRINTER[i+3*nerrors]=mu[jdhk] <meps ? 0  : fabs(mu[jdhk]-1)   <meps ? 1:mu[jdhk];
			TRINTER[i+4*nerrors]=alpha[jdhk]<meps ? 0: fabs(alpha[jdhk]-1)<meps ? 1:alpha[jdhk];
			i++;
		}
	}
	mxFree(dhki);
	mxFree(mu);
	mxFree(alpha);
	mxFree(lamb);
	mxFree(orient);
	mxFree(dhk);
	mxFree(pnt);	
}

/**************************** End of distgraph.c *****************************/


