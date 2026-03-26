/*************************** graphdist.c ******************************************************/
/*                                                                                    
/* GRAPHDIST \nCalculates the shortest paths in a fully connected graph		
/*
/* SYNTAX
/* ADJ=graphdist(VER,ITRI,mode)
/* [ADJ,DIST]=graphdist(VER,ITRI,mode)
/* [ADJ,DIST,PATHS]=graphdist(VER,ITRI,mode)
/* [ADJ]=graphdist(ITRI)
/* [ADJ,DIST]=graphdist(ITRI)
/* [ADJ,DIST,PATHS]=graphdist(ITRI)
/* [DIST]=graphdist(ADJ)
/* [DIST,PATHS]=graphdist(ADJ)
/*
/* DESCRIPTION	
/* ADJ = graphdist(ITRI,VER,mode), Returns the adjacency matrix. 
/*							 The adjacency matrix is determined by the mode.
/*                    The following 2 options to calculate the adjacent distances are available: 
/*                      - mode 1         : Adjacency matrix weighted by the eddge distances.
/*                      - mode 2         : Adjacency type of matrix: the distances used being 
/*														 those between individual nodes, provided that they can 
/*													    ''see'' one another in the interior of  the triangulated 
/*														 surface.
/*						- mode 3         : As mode 2, but now the second order neigbors are treated as a neighbor too
/*						- mode 4         : As mode 1, but now the second order neigbors are treated as a neighbor too
/* [ADJ,DIST] = graphdist(ITRI,VER,mode), Returns the adjacency matrix and the distance matrix. 
/*														The distance matrix contains all distance between the 
/*														different vertices.
/*[ADJ,DIST,PATH] = graphdist(ITRI,VER,mode), Returns the adjacency matrix, the distance matrix, 
/*														and the prominent routes matrix. Each column in the 
/*														prominent routes matrix contains the node use for with 
/*														the vertex with the column number as afocus.(see 
/*														van Dam PM, van Oosterom A:Atrial Excitation Assuming 
/*														Uniform Propagation. J Cardiovasc Electrophysiol. 
/*														2003;14:S166-S171)		
/*ADJ            = graphdist(ITRI), Returns the adjacency matrix.
/*[ADJ,DIST]     = graphdist(ITRI), Returns the adjacency and distance matrix.
/*[ADJ,DIST,PATH]= graphdist(ITRI), Returns the adjacency, distance and prominent routes matrix.
/*ADJ            = graphdist(ADJ),  Returns the adjacency matrix.
/*[ADJ,DIST]     = graphdist(ADJ),  Returns the adjacency and distance matrix.
/*[ADJ,DIST,PATH]= graphdist(ADJ),  Returns the adjacency, distance and prominent routes matrix.
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

#define MAXNODES 10000
#define pi 3.14159265358979323846
/************************* Defintions and globals **************************/

int		nrib,npnt,ndhk;
double	maxlength,eps;
int	    (*dhk)[3]		= NULL;
double  (*pnt)			= NULL;
double  (*ADJ)			= NULL;
double  (*DIST)			= NULL;
double  (*PATH)			= NULL;
char    (*fixed)		= NULL;
int		(*rib)[2]		= NULL;
/***************************************************************************/

void showHelp(void)
{
	mexPrintf("GRAPHDIST Version: 2.02\nCalculates the shortest paths in a fully connected graph\n\n");
	mexPrintf("SYNTAX\n");
	mexPrintf("[DIST,PATHS]=graphdist(ADJ) \n");
	mexPrintf("[ADJ,DIST,PATHS]=graphdist(ITRI) \n");
	mexPrintf("[ADJ,DIST,PATHS]=graphdist(ITRI,VER,mode) \n\n");
	
	mexPrintf("DESCRIPTION\n");		
	mexPrintf("[DIST,PATH]= graphdist(ADJ),  Returns the distance matrix and (optionally) the prominent routes matrix based on \n");
	mexPrintf("                      a precomputed (ADJ) adjacency matrix of a weighted or unweighted graph. \n");
	mexPrintf("                      DIST iS the distance matrix. The distance matrix contains the distances between all nodes, while traveling along edges.\n");
	mexPrintf("                      PATH is the the prominent routes matrix. Each column (K) in PATH contains the number of times the node\n");
	mexPrintf("                      is included in all possible shortest paths starting at node K.(see van Dam PM, van Oosterom A\n");
	mexPrintf("                      Atrial Excitation Assuming Uniform Propagation. J Cardiovasc Electrophysiol. 2003;14:S166-S171)\n\n");		
	mexPrintf("[ADJ,DIST,PATH]= graphdist(ITRI), Returns the (UN-WEIGHTED) adjacency matrix derived pertaining to the triangle specification ITRI \n");
	mexPrintf("                      of a triangulated surface. Optionally the corresponding distance and prominent routes matrix is returned.\n");
	mexPrintf("[ADJ,DIST,PATH]= graphdist(ITRI,VER,mode), As previous one, but now based on a WEIGHTED graph. \n");
	mexPrintf("                     The following 2 options to calculate the involved distances are available: \n");
	mexPrintf("                      - mode 1 (default) : based on edge lengths. \n");
	mexPrintf("                      - mode 2           : based on internodal distantances in 3D space, provided \n");
	mexPrintf("                                           that they can 'see' one another in the interior of the triangulated surface. \n");
	mexPrintf("                      - mode 3           : As mode 2, but now the second order neigbors are treated as direct neighbors too.\n");
	mexPrintf("                      - mode 4           : As mode 1, but now the second order neigbors are treated as direct neighbors too.\n");	
	mexPrintf("See also\n- getroute, graphdistone \n");
}
/***************************************************************************/
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

void veccross(double *r1, double *r2, double *r3)
{
  r3[0]=r1[1]*r2[2]-r1[2]*r2[1];
  r3[1]=r1[2]*r2[0]-r1[0]*r2[2];
  r3[2]=r1[0]*r2[1]-r1[1]*r2[0];
 }

static double deter(double r1[3], double r2[3], double r3[3])
 {  
  return r1[0]*(r2[1]*r3[2]-r2[2]*r3[1])
		-r1[1]*(r2[0]*r3[2]-r2[2]*r3[0])
		+r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]);
}


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
double calcSecondOrderLength(int i1,int i2,int i3, int j2, double alpha)
{
	double v1[3],v2[3],vm;
	int j;
	
	for (j=0;j<3;j++)
	{
		vm=pnt[i2+j*npnt]+(pnt[i3+j*npnt]-pnt[i2+j*npnt])*(alpha);
		v1[j]=vm-pnt[i1+j*npnt];
		v2[j]=vm-pnt[j2+j*npnt];				
	}
	return norm(v1)+norm(v2);
}


double calcLength(int i1,int j2)
{
	int i,j,i2,i3,k;
	double alpha,Dleft,Dright,xleft,xright,Dm,delta=1e-6,dist=-1;

	for (i=0;i<ndhk;i++)
	{
		j=ndhk;
		if (i1==dhk[i][0])
		{	i2=dhk[i][1];i3=dhk[i][2];}
		else if (i1==dhk[i][1] )
		{	i2=dhk[i][2];i3=dhk[i][0];}
		else if (i1==dhk[i][2] )
		{	i2=dhk[i][0];i3=dhk[i][1];}
		else
			continue;

		for (j=0;j<ndhk;j++)
		{
			if		(j2==dhk[j][0] && (i2==dhk[j][1] && i3==dhk[j][2] || i2==dhk[j][2] && i3==dhk[j][1]))
				break;
			else if (j2==dhk[j][1] && (i2==dhk[j][0] && i3==dhk[j][2] || i2==dhk[j][2] && i3==dhk[j][0]))
				break;
			else if (j2==dhk[j][2] && (i2==dhk[j][0] && i3==dhk[j][1] || i2==dhk[j][1] && i3==dhk[j][0]))
				break;
		}
		if (j!=ndhk) /*found both triangles */ 
		{
			xleft=0;xright=1;
			Dleft=ADJ[i1+npnt*i2]+ADJ[j2+npnt*i2];
			Dright=ADJ[i1+npnt*i3]+ADJ[j2+npnt*i3];
			for (k=0;k<20 ;k++)
			{
				alpha=(xleft+xright)/2;
				Dm=calcSecondOrderLength(i1,i2,i3, j2, alpha);
				if ((Dleft-Dm)>(Dright-Dm) )
				{
					xleft=alpha;
					Dleft=Dm;
				}
				else
				{
					xright=alpha;
					Dright=Dm;
				}
			}
			if (dist<0)
				dist=Dm;
			else if (Dm<dist)
				dist=Dm;
		}
	}
	return 	dist;
}

/***************************************************************************/
int nodeview(int j1,int j2)
{
	/* view =1 if line connecting node j1 to node j2 of
	 the triangulated surface specified by (VER and ITRI) lies entirely in 
	 its interior; else view=0*/

	double a1[3],a2[3],a3[3],sav[3],lambda,mu;
	double r1[3],r2[3],r3[3],det,position[3],hoek;
	int i,j,na,ns,p; 
	double alphas[100], deltas[100],alpha,tmin;

	for (i=0;i<3;i++)
		a3[i]=pnt[j1+i*npnt]-pnt[j2+i*npnt];

	ns=0;alphas[0]=0;
	/* identify intersections in between node1 and node2*/
	for (i=0;i<ndhk;i++)
	{
        if ((dhk[i][0]==j1 || dhk[i][1]==j1 || dhk[i][2]==j1) &&
            (dhk[i][0]==j2 || dhk[i][1]==j2 || dhk[i][2]==j2))
            return 1;
        if ((dhk[i][0]==j1 || dhk[i][1]==j1 || dhk[i][2]==j1) ||
            (dhk[i][0]==j2 || dhk[i][1]==j2 || dhk[i][2]==j2))
            continue;
			for (p=0;p<3;p++)
			{
				a1[p]=pnt[dhk[i][1]+p*npnt]-pnt[dhk[i][0]+p*npnt];
				a2[p]=pnt[dhk[i][2]+p*npnt]-pnt[dhk[i][0]+p*npnt];
 			}
			det=deter(a1, a2, a3);
	  	 	if (fabs(det) > eps)
			{
	     	   	/* no paralel detected; solution via Cramer's rule*/
			 	for (p=0;p<3;p++)
					sav[p]=pnt[j1+p*npnt]-pnt[dhk[i][0]+p*npnt];
				lambda=deter(sav,a2,a3)/det;
	 	       	if (lambda >= 0 && lambda <= 1)
		      	{
			    	mu=deter(a1, sav, a3)/det;
					if (mu >= 0 && mu <= 1 && (lambda + mu) <= 1)
				  	{
    		  	      	alpha=deter(a1,a2,sav)/det; 
               	       	if (alpha > 0 && alpha < 1)
						{
							ns++;
                          	alphas[ns]=alpha;            
						}
					}	
				}
		}
	}
	alphas[++ns]=1;	ns++;
	deltas[0]=0;	na=1;
	/* sort for increasing values of alpha*/
	for (i=0;i<ns;i++)
	{
		tmin=1;
		for (j=0;j<ns;j++)
			if (alphas[j]>deltas[na-1] && alphas[j] <tmin)
				tmin=alphas[j];
		if (tmin==1) 
			break;
		deltas[na]=tmin;
		na++;
	}
	deltas[na]=1.0;
	na++;
	
	/* test if positions halfway the intersections are interior or exterior */
	/* if any of these is exterior, nodes j1 and j2 cannot be connected     */
	for (i=1; i<na;i++)
	{
		for(p=0;p<3;p++)
			position[p]=pnt[j1+p*npnt]+((deltas[i]+deltas[i-1])/2)*(pnt[j2+p*npnt]-pnt[j1+p*npnt]);
		
		hoek=0;	
		for (j=0;j<ndhk;j++)
		{
			for(p=0;p<3;p++)
			{
				r1[p]=pnt[dhk[j][0]+npnt*p]-position[p];				
				r2[p]=pnt[dhk[j][1]+npnt*p]-position[p];
				r3[p]=pnt[dhk[j][2]+npnt*p]-position[p];			
			}	
			hoek+=rhoek(r1,r2,r3);
		}
		if (fabs(hoek) < 0.1)
			return 0; 
	}
	return 1;
}

/***************************************************************************/
void wallDist(void)
{
	int j1,j2,i;
	double r[3];
 	for (j1=0;j1<npnt;j1++)
		for (j2=j1+1; j2<npnt;j2++)
		{
			if (nodeview(j1,j2)==1)
			{
				for (i=0;i<3;i++)
					r[i]=pnt[j1+i*npnt]-pnt[j2+i*npnt];
				ADJ[j1+npnt*j2]=norm(r);
				ADJ[j2+npnt*j1]=ADJ[j1+npnt*j2];
			}
		}
}
/***************************************************************************/
void createAdjacencyItri(void)
{
	int i,j,j1,j2;

	for (j1=0;j1<npnt;j1++)
		for (j2=j1+1; j2<npnt;j2++)
		{
			ADJ[j1+npnt*j2]=0;
			ADJ[j2+npnt*j1]=0;
		}	
	
	for (i=0;i<ndhk;i++)
		for (j=0;j<3;j++)
		{
			j1=icyc(j,3);j2=icyc(j+1,3);
			ADJ[dhk[i][j1]+dhk[i][j2]*npnt]=1;
			ADJ[dhk[i][j1]*npnt+dhk[i][j2]]=1;
		}
}

/***************************************************************************/
void createAdjacency2D()
{
	int i,j,j1,j2,t;
	double r[3];

	for (j1=0;j1<npnt;j1++)
		for (j2=j1+1; j2<npnt;j2++)
		{
			ADJ[j1+npnt*j2]=0;
			ADJ[j2+npnt*j1]=0;
		}	
	
	for (t=0;t<ndhk;t++)
	{
		for (j=0;j<3;j++)
		{
			j1=icyc(j,3);j2=icyc(j+1,3);
			for (i=0;i<3;i++)
				r[i]=pnt[dhk[t][j1]+i*npnt]-pnt[dhk[t][j2]+i*npnt];
			ADJ[dhk[t][j1]+dhk[t][j2]*npnt]=norm(r);
			ADJ[dhk[t][j1]*npnt+dhk[t][j2]]=ADJ[dhk[t][j1]+dhk[t][j2]*npnt];
		}
	}
}
/***************************************************************************/
void calcAdjacency3D(void)
{		
	createAdjacency2D();
	wallDist();
}

/***************************************************************************/
void calcAdjacency2DsecondOrder(void)
{
	int j1,i,j,k;
	double rn1[3],rn2[3],v[3],beta,len;
	char *DD=NULL;

	createAdjacency2D();
	DD=mxCalloc(npnt*npnt, sizeof(char));	DD|| outofmem();
	for (i=0;i<npnt;i++)
		for (j=0;j<npnt;j++)
            if (ADJ[i+j*npnt]>0)
    			DD[i+j*npnt]=1;            
            else
    			DD[i+j*npnt]=0;
	for (i=0;i<npnt;i++)
		for (j=0;j<npnt;j++)
			if (DD[i+j*npnt])
			{
				for (k=0;k<npnt;k++)
					if (DD[k+j*npnt] && !DD[i+k*npnt] && k!=i)
					{
						len=calcLength(i,k);
						if (len>0 & (len < ADJ[i+npnt*k] | ADJ[i+npnt*k]==0) )
						{
							ADJ[i+npnt*k]=len;
							ADJ[k+npnt*i]=len;
						}
					}
			}
	mxFree(DD);
}

/*******************************************************************/
void calcAdjacency3DsecondOrder(void)
{
	calcAdjacency2DsecondOrder();
	wallDist();	
}

/***************************************************************************/
void init_getdistmat(void)
{
	int i,j,k;
 	nrib=0;
	for ( i=0;i<npnt;i++)
	 	for (j=i;j<npnt;j++)
		{
			DIST[i*npnt+j]=ADJ[i*npnt+j];
			DIST[j*npnt+i]=ADJ[i*npnt+j];
			if (DIST[i*npnt+j]>0) 
				nrib++;
		}
	rib   = mxCalloc(nrib,2*sizeof(int));   rib	  || outofmem();
	k=0;  maxlength = 0;
	/* create a list of all the ribs */
	for ( i=0;i<npnt;i++)
	 	for (j=i;j<npnt;j++)
    		if (DIST[i*npnt+j]>0)
			{
	 		 	rib[k][0]=i;
				rib[k][1]=j;
				maxlength = maxlength+DIST[i*npnt+j];
				k++;
			}
	fixed = mxCalloc(npnt,sizeof(char));    fixed || outofmem();
	for (i=0; i<npnt;i++)
	{	
		for (j=0;j<npnt;j++)
		{
			if (DIST[i+j*npnt]==0 )
				DIST[i+j*npnt]=maxlength;
		}
		fixed[i]=0;
	}
}

/***************************************************************************/
void updateRibs(int start)
{
	int i,k=0;

	for ( i=0;i<nrib;i++)
	{
		if (rib[i][0]>=start || rib[i][1]>=start) 
		{
			rib[k][0]=rib[i][0];
			rib[k][1]=rib[i][1];		
			k++;
		}
	}	
	nrib=k;
}
/***************************************************************************/
void CalculateDists(void)
{
	double length;
	int fixrib[MAXNODES],i,n,nfix,start, maxNodes=MAXNODES-1;

	init_getdistmat();			
	for (start=0;start<npnt;start++)
	{
		nfix=0;
		updateRibs(start);		
		DIST[start+start*npnt]=0;
		for (i=0; i<npnt;i++)
			if (i<=start)
			{
				fixed[i]=1;
				nfix++;
			}
			else
				fixed[i]=0;

		/* calculate all paths between startpnt and vertices*/
		while (nfix<npnt)
		{
			length = maxlength-0.01;
			n=0;
			for (i=0;i<nrib;i++)
			{
				if ( !fixed[rib[i][1]] && fixed[rib[i][0]] &&
					DIST[rib[i][0]+start*npnt]+ DIST[rib[i][0]+rib[i][1]*npnt]<=length)
				{
					length = DIST[rib[i][0]+start*npnt]+ DIST[rib[i][0]+rib[i][1]*npnt];
					fixrib[n]=i;
					n++;
				}
				else if(!fixed[rib[i][0]] && fixed[rib[i][1]] &&
						DIST[rib[i][1]+start*npnt]+ DIST[rib[i][1]+rib[i][0]*npnt] <=length)
				{
					length = DIST[rib[i][1]+start*npnt]+ DIST[rib[i][1]+rib[i][0]*npnt];
					fixrib[n]=i;
					n++;
				}
			}
			if (n==0)
				mexErrMsgTxt("Inconsistent graph (more than 1 part).");
			if (n>maxNodes)
				mexErrMsgTxt("Too many nodes (>1000) with the same length.");
			for (i=0;i<n;i++)
			{
				if( fixed[rib[fixrib[i]][1]] &&
					 fabs(DIST[rib[fixrib[i]][1]+start*npnt]+ DIST[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < eps)      
				{
					fixed[rib[fixrib[i]][0]]=1;
					DIST[rib[fixrib[i]][0]+start*npnt]=length;
					DIST[start+rib[fixrib[i]][0]*npnt]=length;
				}
				else if (fabs(DIST[rib[fixrib[i]][0]+start*npnt]+ DIST[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < eps)      
				{
					fixed[rib[fixrib[i]][1]]=1;
					DIST[rib[fixrib[i]][1]+start*npnt]=length;
					DIST[start+rib[fixrib[i]][1]*npnt]=length;
				}
			}
			nfix=0;
			for (i=0;i<npnt;i++)
            {
				if (fixed[i]) 
					nfix++;
            }
		}
	}
	mxFree(fixed);
	mxFree(rib);
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
		DIST[start+start*npnt]=0;
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
						DIST[rib[i][0]+start*npnt]+ DIST[rib[i][0]+rib[i][1]*npnt]<=length)
				{
					length = DIST[rib[i][0]+start*npnt]+ DIST[rib[i][0]+rib[i][1]*npnt];
					fixrib[n]=i;
					n++;
				}
				else if(!fixed[rib[i][0]] && fixed[rib[i][1]] &&
							DIST[rib[i][1]+start*npnt]+ DIST[rib[i][1]+rib[i][0]*npnt] <=length)
				{
					length = DIST[rib[i][1]+start*npnt]+ DIST[rib[i][1]+rib[i][0]*npnt];
					fixrib[n]=i;
					n++;
				}
			}
			if (n==0)
				mexErrMsgTxt("Inconsistent graph (more than 1 part).");
			if (n>maxNodes)
				mexErrMsgTxt("Too many nodes (>1000) with the same length.");
			for (i=0;i<n;i++)
			{
				if(fixed[rib[fixrib[i]][1]] &&
					 fabs(DIST[rib[fixrib[i]][1]+start*npnt]+ DIST[rib[fixrib[i]][1]+rib[fixrib[i]][0]*npnt]- length) < eps)      
				{
					fixed[rib[fixrib[i]][0]]=1;
					DIST[rib[fixrib[i]][0]+start*npnt]=length;
					DIST[start+rib[fixrib[i]][0]*npnt]=length;
					frompnt[rib[fixrib[i]][0]]=rib[fixrib[i]][1];
				}
				else if (fabs(DIST[rib[fixrib[i]][0]+start*npnt]+ DIST[rib[fixrib[i]][0]+rib[fixrib[i]][1]*npnt]- length) < eps)      
				{
					fixed[rib[fixrib[i]][1]]=1;
					DIST[rib[fixrib[i]][1]+start*npnt]=length;
					DIST[start+rib[fixrib[i]][1]*npnt]=length;
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
				PATH[j+start*npnt]++;
				j=frompnt[j];
			}
		}
		PATH[start+start*npnt]=npnt;
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
	int  mrows,mcols,mode,i,j,dodist=0,dopath=0;
	double *itri;
	/* initialize globals*/
	eps=epsilon();
	ndhk=0;
	npnt=0;
	nrib=0;
	/* Check input arguments */
	if(nrhs ==2 || nrhs >3 || nlhs >3 ) 
	{
		showHelp();
		mexErrMsgTxt("Inproper number of arguments ");
	}
	if (nrhs ==3 || nrhs==2) /* ITRI VER (mode) input */
	{
		if (nlhs>1) dodist=1;
		if (nlhs>2) dopath=1;
		mode=1;		/* default only surface */
		/* ITRI */ 
		mrows=mxGetM(prhs[0]);
	 	mcols=mxGetN(prhs[0]);  	 	
		if( mcols !=3)				mexErrMsgTxt("size ITRI = [.., 3]");
		ndhk = mrows;
		itri = mxGetPr(prhs[0]);
		for (i=0;i<ndhk;i++)
			if (itri[i]-(int)itri[i]!=0 )
				mexErrMsgTxt("first parameter (ITRI) must contain index numbers only");
				
		/* VER	*/
	 	mrows=mxGetM(prhs[1]);
		mcols=mxGetN(prhs[1]);  	 	
		if( mcols !=3)				mexErrMsgTxt("size VER = [.., 3]");
		if(!(mxIsDouble(prhs[1]))) 	mexErrMsgTxt("VER array must be of type double.");
		npnt = mrows;
		pnt  = mxGetPr(prhs[1]);

		/* mode */
		if (nrhs==3)
		{
			if( mxGetM(prhs[2])!=1 && mxGetN(prhs[2])!=1 )	mexErrMsgTxt("mode must be an scalar");
			mode=(int)mxGetScalar(prhs[2]);
		}
		if (mode>4 || mode<1)
		{
			showHelp();
			mexErrMsgTxt("mode must be between 1 and 4");
		}
		dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();
		for (i=0;i<ndhk;i++)
		{
			dhk[i][0]=(int)itri[i]-1;
			dhk[i][1]=(int)itri[i+ndhk]-1;		
			dhk[i][2]=(int)itri[i+2*ndhk]-1;	
			if (dhk[i][0]>npnt ||dhk[i][1]>npnt || dhk[i][2]>npnt)
			{
				mexErrMsgTxt("Index numbers in ITRI are larger than number of vertices in VER");
				mxFree(dhk);
			}
		}
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);					ADJ  = mxGetPr(plhs[0]);
		if (dodist)	{ plhs[1]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);	DIST = mxGetPr(plhs[1]);}
		if (dopath)	{ plhs[2] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);	PATH = mxGetPr(plhs[2]);}		
		if (mode==1 )
			createAdjacency2D();
		else if (mode==2)
			calcAdjacency3D();	
		else if (mode==3)
			calcAdjacency3DsecondOrder();
		else if (mode==4)
			calcAdjacency2DsecondOrder();
		mxFree(dhk);
	}
	else if (nrhs ==1 && mxGetM(prhs[0])!= mxGetN(prhs[0])) //not a square matrix so ITRI
	{
		if (nlhs>1) dodist=1;
		if (nlhs>2) dopath=1;

		/* ITRI */
	 	mrows=mxGetM(prhs[0]);	ndhk = mrows;
		mcols=mxGetN(prhs[0]);  if( mcols !=3)	mexErrMsgTxt("size ITRI = [.., 3]");
		itri = mxGetPr(prhs[0]); 
		dhk=mxCalloc(ndhk, 3*sizeof(int));	dhk	|| outofmem();
		for (i=0;i<ndhk;i++)
		{
			dhk[i][0]=(int)itri[i]-1;
			dhk[i][1]=(int)itri[i+ndhk]-1;		
			dhk[i][2]=(int)itri[i+2*ndhk]-1;	
			if (dhk[i][0]>npnt) npnt=dhk[i][0];
			if (dhk[i][1]>npnt) npnt=dhk[i][1];
			if (dhk[i][2]>npnt) npnt=dhk[i][2];
		}
		npnt++;
		if (npnt < ndhk/5)
		{
			mxFree(dhk);
			mexErrMsgTxt("number vertices is less than 1/5 of the number of triangles, Possible wrong input? ");
		}
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);	ADJ = mxGetPr(plhs[0]);
		
		createAdjacencyItri();
		if (dodist)	{plhs[1] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);DIST = mxGetPr(plhs[1]);}
		if (dopath)	{plhs[2] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);PATH = mxGetPr(plhs[2]);}		
		mxFree(dhk);
	}
	else if (nrhs==1) /* adjacency matrix input */
	{
		dodist=1;
		if (nlhs>1) dopath=1;
	 	mrows=mxGetM(prhs[0]);		mcols=mxGetN(prhs[0]);  		npnt = mrows;
		if( mcols !=mrows)			mexErrMsgTxt("The adjacency matrix must be square");
		if(!(mxIsDouble(prhs[0]))) 	mexErrMsgTxt("ADJ array must be of type double.");
		ADJ  = mxGetPr(prhs[0]);
		
		/*Set output pointer to output matrix*/   
		plhs[0]=mxCreateDoubleMatrix(npnt,npnt,mxREAL);					DIST = mxGetPr(plhs[0]);
		if (dopath){plhs[1] = mxCreateDoubleMatrix(npnt,npnt,mxREAL);	PATH = mxGetPr(plhs[1]);}
	}
	else
		showHelp();
	
	
	if (dodist)
		if ( dopath )
			CalculateDistsAndPaths();
		else 
			CalculateDists();
}

/**************************** End of graphdist.c *****************************/
