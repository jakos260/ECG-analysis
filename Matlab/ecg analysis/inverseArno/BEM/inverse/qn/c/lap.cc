/***************************************************************************/
/* Program lap                                                             */
/* Geertjan Huiskamp                                                       */
/* Medical physics and biophysics                                          */
/* University of Nijmegen The Netherlands, UCI june '92                    */
/*                                  t                                      */
/* Calculate regularisation-matrix C C for Laplacian operator              */
/* Laplacian estimator with angular weightfactors                          */
/* Area normalization                                                      */
/* Functions used from files: matrix.c, nrutils.c,                         */
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trilib.h"
#define NNBMAX 257

/*
double *dvector(int,int);
double **dmatrix(int,int,int,int);
void   free_dvector(double *,int,int);
void   free_dmatrix(double **,int,int,int,int);

extern double **matrix(int,int,int,int);
extern int   *ivector(int,int);
extern int   **imatrix(int,int,int,int);
extern void  ifree_vector(int *,int,int);
extern void  ifree_matrix(int **,int,int,int,int);
extern int   binputmatrix(char *,double **,int,int);
*/

extern int voisin(int **,int,int,int ***,int **,int **,int ***,int **,int);

/*
int dhkinp(char *filename,int *npnt,int *ndhk,double (**pnt1)[3],int (**dhk1)[3]);
*/

double *neigharea(int,int,double **,int **);
void  distanc(int,int *,int **,double **,double *,double **);
void  angwei(int,int *,int **,double **,double **,double *,double **);
double **cmat(int,int *,int **,double *,double *,double *,double **,double **,int);
double **ctcmat(int,int *,int **,double **);
void  outctc(char *,int,double **);

int main (int argc, char *argv[])      
{
 int   i,j,k,idum,npth,ntrh,maxnnb;
 int   **itrh,*nnb,**nb,*nnbtr,**nbtr,*opn,(*dhk)[3];
 double **pnth,**doublemat,(*pnt)[3];
 double *area,*dm,*wtot,**d,**wgt,**c,**ctc;
 char  heartname[80],outname[80];
 FILE  *hfp;

 setbuf(stdout,NULL);

 printf("triagulated surface>");
 gets(heartname);

/*
 hfp=fopen(heartname,"r");
 fscanf(hfp,"%d\n",&npth);
 pnth=matrix(1,npth,1,3); 
 for(i=1;i<=npth;i++) 
 fscanf(hfp,"%d %lf %lf %lf\n",&idum,&pnth[i][1],&pnth[i][2],&pnth[i][3]);
 fscanf(hfp,"%d\n",&ntrh);
 itrh=imatrix(1,ntrh,1,3);
 for(i=1;i<=ntrh;i++) 
 fscanf(hfp,"%d %d %d %d\n",&idum,&itrh[i][1],&itrh[i][2],&itrh[i][3]);
 close(hfp);
*/

  if (!dhkinp(heartname, &npth, &ntrh, &pnt, &dhk))
    exit(1);
 pnth=matrix(1,npth,1,3); 
 for(i=1;i<=npth;i++)
   for (k=1; k<=3; k++)
     pnth[i][k]=pnt[i-1][k-1];
 itrh=imatrix(1,ntrh,1,3);
 for(i=1;i<=ntrh;i++) 
   for (k=1; k<=3; k++)
     itrh[i][k]=dhk[i-1][k-1]+1;

 maxnnb=voisin(itrh,ntrh,npth,&nb,&nnb,&opn,&nbtr,&nnbtr,NNBMAX);
 printf("maximum number of neighbours %d\n",maxnnb);

 area=neigharea(npth,ntrh,pnth,itrh);

 ifree_matrix(itrh,1,ntrh,1,3);

 dm=vector(1,npth);
 d=matrix(1,npth,1,maxnnb);
 distanc(npth,nnb,nb,pnth,dm,d);

 wtot=vector(1,npth);
 wgt=matrix(1,npth,1,maxnnb);
 angwei(npth,nnb,nb,pnth,d,wtot,wgt);
 free_matrix(pnth,1,npth,1,3);

 c=cmat(npth,nnb,nb,area,wtot,dm,wgt,d,1);
 printf("\noutput Laplacian matrix>");
 gets(outname);

 if(outname[0]!='\0') outctc(outname,npth,c);

 free_matrix(d,1,npth,1,maxnnb);
 free_matrix(wgt,1,npth,1,maxnnb);

 printf("\noutput LtL matrix>");
 gets(outname);
 if(outname[0]=='\0') exit(1);

 ctc=ctcmat(npth,nnb,nb,c);
 free_matrix(c,1,npth,1,npth);
 ifree_vector(nnb,1,npth);
 ifree_matrix(nb,1,npth,0,NNBMAX);

 outctc(outname,npth,ctc);
 free_matrix(ctc,1,npth,1,npth);
 return 0;
}
 

double *neigharea(int npth,int ntrh,double **pnth,int **itrh)
{
 int     i,j,i1,i2,i3;
 double  *area,rm[4],rp[4],ov1,ov2,ov3,opp;

 area=vector(1,npth); 
 for(i=1;i<=npth;i++) area[i]=0.0;
 for(i=1;i<=ntrh;i++)
    {
     i1=itrh[i][1];
     i2=itrh[i][2];
     i3=itrh[i][3];
     for(j=1;j<=3;j++)
        {
         rm[j]=pnth[i2][j]-pnth[i1][j];
         rp[j]=pnth[i3][j]-pnth[i1][j];
        }
     ov1=rp[2]*rm[3]-rm[2]*rp[3];
     ov2=rp[3]*rm[1]-rm[3]*rp[1];
     ov3=rp[1]*rm[2]-rm[1]*rp[2];
     opp=sqrt(ov1*ov1+ov2*ov2+ov3*ov3)/6;
     area[i1]+=opp;
     area[i2]+=opp;
     area[i3]+=opp;
    }
 return area;
}

void distanc(int npth,int *nnb,int **nb,double **pnth,double *dm,double **d)
{
 int    i,j,k;
 double  tmp;

 for(i=1;i<=npth;i++)
    {
     dm[i]=0;
     for(j=1;j<=nnb[i];j++)
        {
         tmp=0;
         for(k=1;k<=3;k++)
         tmp+=(pnth[nb[i][j]][k]-pnth[i][k])*(pnth[nb[i][j]][k]-pnth[i][k]);
         d[i][j]=sqrt(tmp);
         dm[i]+=d[i][j];
        }
      dm[i]=dm[i]/nnb[i];
    }
}

void angwei(int npth,int *nnb,int **nb,double **pnth,double **d,double *wtot,double **wgt)
{
 int   i,j,k,jm,jp;
 double ro[4],rm[4],rp[4],argm,argp,apls,amin,pim,pip,tmp;

 for(i=1;i<=npth;i++)
    {
     wtot[i]=0;
     for(j=1;j<=nnb[i];j++)
        {
         jm=j-1;
         if(jm==0) jm=nnb[i];
         jp=j+1;
         if(jp>nnb[i]) jp=1;
         pim=0;
         pip=0;
         for(k=1;k<=3;k++)
            {
             ro[k]=pnth[nb[i][j]][k]-pnth[i][k];
             rm[k]=pnth[nb[i][jm]][k]-pnth[i][k];
             rp[k]=pnth[nb[i][jp]][k]-pnth[i][k];
             pim+=rm[k]*ro[k];
             pip+=rp[k]*ro[k];
            }
         argm=pim/(d[i][jm]*d[i][j]);
         if(argm>1) argm=1.0; 
         else if(argm<-1) argm=-1.0; 
         amin=acos(argm);
         argp=pip/(d[i][j]*d[i][jp]);
         if(argp>1) argp=1.0;
         else if(argp<-1) argp=-1.0; 
         apls=acos(argp);
         wgt[i][j]=tmp=((1-argm)/sin(amin))+((1-argp)/sin(apls));
         wtot[i]+=tmp;
        }
    }
}

double **cmat(int npth,int *nnb,int **nb,double *area,double *wtot,double *dm,double **wgt,double **d,int warea)
{
 int   i,j;
 double **c,tmp,fac; 

 printf("setup Laplacian matrix"); 
 if(warea) printf(", area weighted\n"); 
 c=matrix(1,npth,1,npth);
 for(i=1;i<=npth;i++) for(j=1;j<npth;j++) c[i][j]=0.0;
 for(i=1;i<=npth;i++)
    {
     if(warea) fac=4*sqrt(area[i])/wtot[i]/dm[i]; else fac=4/wtot[i]/dm[i];
     tmp=0;
     for(j=1;j<=nnb[i];j++)
        {
         c[i][nb[i][j]]=fac*wgt[i][j]/d[i][j];
         tmp-=c[i][nb[i][j]];
        }
     c[i][i]=tmp;
    }
 return c;
}

double **ctcmat(int npth,int *nnb,int **nb,double **c)
{
 int i,j,k,l,nbik,nbjl;
 double **ctc,tmp,check;
 
 ctc=matrix(1,npth,1,npth);
 printf("     :%d",npth);
 for(i=1;i<=npth;i++)
    {
     printf ("\r %d",i);
     for(j=i;j<=npth;j++)
        {
         tmp=0;
         for(k=0;k<=nnb[i];k++)
            {
             nbik=nb[i][k];              
             for(l=0;l<=nnb[j];l++) 
                { 
                 nbjl=nb[j][l];
                 if(nbik==nbjl) 
		   {
                    tmp+=c[nbik][i]*c[nbjl][j]; 
                    break;                                  
                   }
                }
            }
         ctc[i][j]=tmp;
         ctc[j][i]=tmp;
        }
    } 
 return ctc;
}

void outctc(char *outname,int npth,double **ctc)
{
 int i,j,nnb,maxnnb=0,*index;
 double out;
 FILE *outfile;

 outfile=fopen(outname,"w");
 index=ivector(1,npth);
 fwrite(&npth,sizeof(int),1,outfile); 
 for(i=1;i<=npth;i++)
    {
     nnb=0;
     for(j=1;j<=npth;j++)
        {
         if(ctc[i][j]!=0)
	   {
	    nnb++;	 
            index[nnb]=j;
	   }
        }
     fwrite(&nnb,sizeof(int),1,outfile); 
     if(nnb>maxnnb) maxnnb=nnb;
     for(j=1;j<=nnb;j++) 
        {
         fwrite(&index[j],sizeof(int),1,outfile);         	 
         out=ctc[i][index[j]];
         fwrite(&out,sizeof(double),1,outfile);         	 
        }	 
    } 
 printf("\nmax. # prim.+sec. neighbors %d\n",maxnnb);
 fclose(outfile); 
 ifree_vector(index,1,npth);
}
/*
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
 double **m;
 int i;

 m=(double **) calloc(nrh-nrl+1, sizeof(double));
 if(!m)
   {
    fprintf(stderr,"\n Sorry, out of memory\n\n\a\a");
    return NULL;
   }
 m-= nrl;
 m[nrl]=(double *) calloc(((nrh-nrl+1)*(nch-ncl+1)), sizeof(double));
 if(!m[nrl])
   {
    fprintf(stderr, "\n Sorry, out of memory\n\n\a\a");
    free(m[nrl]);
    free(m+nrl);
    return NULL;
   }
 m[nrl]-= ncl;
 for(i=nrl+1;i<=nrh;i++) m[i]=m[nrl]+(i-nrl)*(nch-ncl+1);
 return m;
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
 int i;

 for(i=nrh;i>=nrl;i--) free((char *)(m[i]+ncl));
 free((char *)(m+nrl));
}

double *dvector(int nl,int nh)
{
 double *v;

 v=(double *) calloc(nh-nl+1,sizeof(double));
 if(v==NULL) 
  {
   fprintf(stderr,"\n Sorry, out of memory\n\n\a\a");
   return NULL;
  }
 return (v-nl);
 }

void free_dvector(double *v,int nl,int nh)
{
 if (v+nl!=NULL) free((char*) (v+nl));
}
*/ 



