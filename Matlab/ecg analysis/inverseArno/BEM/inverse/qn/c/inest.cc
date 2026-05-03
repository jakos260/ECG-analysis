/***************************************************************************/
/* Program in                                                              */
/* Maria Orgonasova Geertjan Huiskamp                                      */
/* Medical physics and biophysics                                          */
/* University of Nijmegen The Netherlands, November 1988 (.c december l991)*/
/*                                                                         */
/* Program to calculate initial estimate for activation times              */
/* SVD/truncated pseudo-inverse of transfer matrix shifting of tau[x]      */
/* values to minimum residu of simulated and measured potentials           */
/* Functions used from files: inest.c, calreg.c, sigdec.c, vcalc.c,        */
/*                            floggt1.c,svdcmp.c, tau.c,                   */
/*                            routines.c,                                  */
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include "trilib.h"

#define BIG 1e+34

#define max(a,b) 		(a<b ? b : a)
#define min(a,b) 		(a>b ? b : a)

/*
extern double *vector(int,int);
extern int   *ivector(int,int);
extern double **getmatrix(char *,int *,int *);
extern double **matrix(short,short,short,short);
extern int   **imatrix(short,short,short,short);
extern void  free_matrix(double **,short,short,short,short);
*/

extern int   neighbor(int,int,int **,int *,int **);
extern void  svdgj(double **,int,int,double *,double **);
extern void  svdcmp(double **,int,int,double *,double **);
extern void  sigdec(int,int,double **,double *,double **);
extern double **vcalcu(int,int,int,int,double,double *,double **);
extern int   locsup(double,int,int *,int **,double *,int *,double *,double *);

double *rimp(int,int,double **);
void  inctc(char *,int ,int *,int **,double **);
void  trpsinv(int,int,int,double **,double *,double **,double **);
void  calctau(int,int,double **,double *,double *,double *,double *);
//int   shift(int,int,double *,double *,double *,int *);
int   shift(int,int,int,double *,double *,double *,int *);
double tauctcatau(int,int *,int **,double **,double *);
int   optims(int,int,int,int,double **,double **,double,double *);
double referr(int,int,int,double *,double *);

int dhkinp(char *filename,int *npnt,int *ndhk,double (**pnt1)[3],int (**dhk1)[3]);

int main(int argc, char *argv[])
{
 int i,j,k,npntsh,ntrsh,idum,maxnnb,nel,iwind,itref,mxm,iwidth,irank,ibeg,iend;
 int nat,ishift,nsup1,nsup,noref,nnbmax=50, (*dhk)[3];
 int **itrh,*nnb,**nb,*ncontrib,**index,*supind;
 double **a,**u,**v,*sig,**vm,*vin,**ctca,*tau,*taumin, (*pnt)[3];
 double **psi,*winm,*winp,**vi,*tref,*supval,*supwidth;
 double fdum,tmp,src,str,resn,treg,rest,tres,sigmax,tmin,tmax,norm;
 char heartname[80],amaname[80],bsmname[80],refname[80],ctcname[80];
 char qout[80],number[10],outname[80],line[80];
 FILE *fp,*flog=NULL;

 setbuf(stdout,NULL);

 printf("triangulation heart>");
 gets(heartname);

/*
 fp=fopen(heartname,"r");
 fscanf(fp,"%d\n",&npntsh);
 for(i=1;i<=npntsh;i++) 
 fscanf(fp,"%d %lf %lf %lf\n",&idum,&fdum,&fdum,&fdum);
 fscanf(fp,"%d\n",&ntrsh);
 itrh=imatrix(1,ntrsh,1,3);
 for(i=1;i<=ntrsh;i++) 
 fscanf(fp,"%d %d %d %d\n",&idum,&itrh[i][1],&itrh[i][2],&itrh[i][3]);
 fclose(fp);
*/

  if (!dhkinp(heartname, &npntsh, &ntrsh, &pnt, &dhk, NULL))
    exit(1);
/*
 pnth=matrix(1,npntsh,1,3); 
 for(i=1;i<=npntsh;i++)
   for (k=1; k<=3; k++)
     pnth[i][k]=pnt[i-1][k-1];
*/
 itrh=imatrix(1,ntrsh,1,3);
 for(i=1;i<=ntrsh;i++) 
   for (k=1; k<=3; k++)
     itrh[i][k]=dhk[i-1][k-1]+1;

 nnb=ivector(1,npntsh);
 nb=imatrix(1,npntsh,0,nnbmax);
 maxnnb=neighbor(npntsh,ntrsh,itrh,nnb,nb);
 supind=ivector(1,npntsh);
 supval=vector(1,npntsh);
 supwidth=vector(1,npntsh);

 printf("transfer matrix>");
 gets(amaname);
 a=getmatrix(amaname,&nel,&npntsh);

 printf("Source  >");
 gets(line);
 if(line[0]!='\0') sscanf(line,"%lf",&src); else src=40.0;
 str=src/40.0;
 printf("Window  >");
 gets(line);
 if(line[0]!='\0') sscanf(line,"%d",&iwind); else iwind=11;
 printf("interval>");
 gets(line);
 if(line[0]!='\0') sscanf(line,"%d",&itref); else itref=100;

 tref=vector(1,npntsh);
 noref=1;
 printf("reftimes>");
 gets(refname);
 if(refname[0]!='\0')
   { 
    fp=fopen(refname,"r");
    fscanf(fp,"%d\n",&idum);
    if(idum!=npntsh) {printf("\nerror in file %s\n",refname); exit(1);} 
    for(i=1;i<=npntsh;i++) fscanf(fp,"%lf\n",&tref[i]);
    fclose(fp);
    noref=0;
   }

 printf("measurement matrix>");
 gets(bsmname);
 vm=getmatrix(bsmname,&nel,&nat);
 resn=.0;
 for(i=1;i<=nel;i++) for(j=1;j<=itref;j++) resn+=vm[i][j]*vm[i][j];
 vin=rimp(nel,itref,vm);

 printf("regularization matrix>");
 gets(ctcname);
 index=imatrix(1,npntsh,1,nnbmax);
 ncontrib=ivector(1,npntsh);
 ctca=matrix(1,npntsh,1,nnbmax);
 inctc(ctcname,npntsh,ncontrib,index,ctca);

 mxm=max(nel,npntsh);
 u=matrix(1,nel,1,npntsh); 
 v=matrix(1,npntsh,1,npntsh);
 sig=vector(1,mxm);
 for(i=1;i<=nel;i++) for(j=1;j<=npntsh;j++) 
   {
    a[i][j]*=str;
    u[i][j]=a[i][j];
   }
 printf ("Singular value decomposition\n");
 svdgj(u,nel,npntsh,sig,v); 
 sigdec(nel,npntsh,u,sig,v); 

 sigmax=sig[1];


 tau=vector(1,npntsh);
 taumin=vector(1,npntsh);
 winm=vector(1,npntsh);
 winp=vector(1,npntsh);
 psi=matrix(1,npntsh,1,nel);
 
 ibeg=1;
 iend=min(nel,npntsh);

 irank=ibeg;
 
mainloop:

iwidth=1;
if (ibeg!=iend)
  { 
   printf("\noutput monitor file>");
   gets(qout);
   flog=fopen(qout,"w");
   fprintf(flog,"%s initial estimate logfile\n",amaname);
   fprintf(flog,"window %d source %lf norm %d reference %s\n",iwind,src,itref,refname);
   fprintf(flog,"ran\ts/smax\twidth\treg\tres\ttres\tdesc\tnsup\n");
   printf ("\nran  smax  tmin  tmax width      reg     res    tres    desc nsup\n");
  }

 while (irank<=iend && iwidth-2*((iwind-1)/2)<itref && sig[irank]!=0)	 
  {
   trpsinv(irank,nel,npntsh,u,sig,v,psi);            
   calctau(nel,npntsh,psi,vin,tau,&tmin,&tmax);
   norm=0; 
   for(i=1;i<=npntsh;i++) norm+=tau[i]*tau[i];
   fdum=0;
   for(i=1;i<=nel;i++) 
      {
       tmp=0;
       for(j=1;j<=npntsh;j++) tmp+=tau[j]*a[i][j];
       tmp-=vin[i];
       fdum+=tmp*tmp;
      }
   fdum=sqrt(fdum/nel);
   //   if(shift(npntsh,itref,tau,&tmin,&tmax,&iwidth))
   if(shift(npntsh,itref,iwind,tau,&tmin,&tmax,&iwidth))
     {
      treg=tauctcatau(npntsh,ncontrib,index,ctca,tau);
      vi=vcalcu(nel,npntsh,iwidth,iwind,str,tau,a);    
      ishift=optims(nel,itref,iwidth,iwind,vm,vi,resn,&rest);
      free_matrix(vi,1,nel,1,iwidth);  
      tres=referr(noref,npntsh,ishift,tau,tref);     
      nsup=locsup(1.0,npntsh,nnb,nb,tau,supind,supval,supwidth);
      if(ibeg!=iend)
        {
          fprintf(flog,"%3d\t%6.4f\t%6.4f\t%8.2f\t%6.4f\t%8.4f\t%8.4f\t%d\n",irank,
                       sig[irank]/sigmax,norm,treg,rest,tres,fdum,nsup);
          printf("%3d %5.3f %5.2f %5.2f   %3d%10.2f%8.4f%8.4f%8.4f  %2d\n",irank,
	          sig[irank]/sigmax,tmin+ishift,tmax+ishift,(int)(ceil((tmax-tmin))),treg,rest,tres,fdum,nsup);
        }     
      tmin+=ishift;
      tmax+=ishift;
      if(ibeg==iend) 
        {
         printf("output tau(x) file>");
         gets(outname);
         fp=fopen(outname,"w");
         if(fp!=NULL)
           {
            fprintf (fp,"%4d 1\n",npntsh);
            for (i=1;i<=npntsh;i++) fprintf(fp," %lf\n",tau[i]);
            fprintf(fp," case : %s\n",amaname);
            fprintf(fp," #pntsh : %d\n",npntsh);
            fprintf(fp," tmin : %lf\n",tmin);
            fprintf(fp," tmax : %lf\n",tmax);
            fprintf(fp," source: %lf\n",src);
            fprintf(fp," ran/it: %d\n",irank);
            fprintf(fp," reg : %lf\n",treg);
            fprintf(fp," res : %lf\n",rest);
           }
          else
           {
            printf ("Can't open file %s \n",outname);
            exit (1);
           }
          printf ("\nTimes are in file %s \n",outname);  
          fclose (fp);
         }  
        }
      irank++;
  }               

 printf ("\n");

 fclose (flog);
 printf ("Monitoring values are in file %s \n",qout);
 printf("solution at rank:");
 gets(line);
 if(line[0]=='\0') exit(0);
 sscanf(line,"%d",&irank);
 ibeg=irank;
 iend=irank;
 if(irank<=nel && irank >0) goto mainloop;   
 return 0;
}

double *rimp(int nel,int itref,double **vm)
{
 int i,j;
 double *vin; 

 vin=vector(1,nel);
 for(i=1;i<=nel;i++)
    {
     vin[i]=vm[i][1]/2;  
     for(j=2;j<=itref;j++) vin[i]-=(vm[i][j-1]+vm[i][j])/2;
    }
 return vin;
 }

void calctau(int nel,int npntsh,
             double **psi,double *vin,double *tau,double *tmin,double *tmax)
{
 int i,j;

 *tmax=-BIG;
 *tmin=BIG;
 for(i=1;i<=npntsh;i++)
    {
     tau[i]=0;
     for(j=1;j<=nel;j++) tau[i]+=psi[i][j]*vin[j];
     if(tau[i]>*tmax) *tmax=tau[i];
     if(tau[i]<*tmin) *tmin=tau[i];
    }
}

//int shift(int npntsh,int itref,double *tau,double *tmin,double *tmax,int *iwidth)
int shift(int npntsh,int itref,int iwind, double *tau,double *tmin,double *tmax,int *iwidth)
{
 int i,j,left;
 double width;

 left=(iwind-1)/2;
 for(i=1;i<=npntsh;i++) tau[i]-=(*tmin)-1-left;
 *tmax-=(*tmin)-1-left;
 *tmin=1+left;
 width=*tmax-*tmin+2*left;
 *iwidth=((int)(width))+2;
 if(*iwidth<itref+2*left)
   return 1;
 else
   {
    printf("times out of range %lf-%lf\n",*tmin-left,*tmax-left);
    return 0;
   }
}

int optims(int nel,int itref,int iwidth, int iwind,
           double **vm,double **vi,double resn,double *rest)
{
 int is,ie,it,iv,iw,ishift=0,left;
 double resmin,res,vv,ww;

 left=(iwind-1)/2;
 resmin=BIG;
 for (is=1-left;is<=(itref-iwidth+2*left);is++)
 {
  res=0;
  for (ie=1;ie<=nel;ie++)
   for (it=1;it<=itref;it++)
   {
    vv=vm[ie][it];
    iw=it-is+1;
    ww=0;
    if (iw>=0 && iw<=iwidth)
       ww=vi[ie][iw];
    res+=(vv-ww)*(vv-ww);
    }
 if (res<resmin)
 {
  resmin=res;
  ishift=is-1;
  }
 }
 *rest=sqrt((double)(resmin/resn));
 return ishift;
}

double referr(int noref,int npntsh,int ishift,double *tau,double *tref)
{
 int i;
 double tres,trefn,errmax,tmp;

 tres=0;
 trefn=0;
 errmax=0;
 for(i=1;i<=npntsh;i++)
    {
     tau[i]+=ishift;
     trefn+=tref[i]*tref[i];
     tmp=fabs(tref[i]-tau[i]);
     tres+=tmp*tmp;
     if(tmp>errmax) errmax=tmp;
     if(noref) tref[i]=tau[i];
    }
 if(trefn==0) trefn=1;
 tres=sqrt((double)(tres/trefn));
 return tres;
}

void trpsinv(int irank,int nrow,int ncol,
            double **u,double *sig,double **v,double **psi)
{
 int i,j,k,nsig;
 double dum;
 static int lastrank=1;

 nsig=min(nrow,ncol);
 if(irank>nsig) {irank=nsig; printf("rank set at %d\n",irank);}
 if(irank<lastrank) 
   {
    lastrank=1; 
    for(i=1;i<=ncol;i++) for(j=1;j<=nrow;j++) psi[i][j]=0;
    printf("trpsinv initialized\n");
   }
 for(i=1;i<=ncol;i++)
    {
     for(j=1;j<=nrow;j++)
        {
         dum=0;
         for(k=lastrank;k<=irank;k++) dum+= v[i][k]*u[j][k]/sig[k];
         psi[i][j]+=dum;
        }
    }
 lastrank=irank+1;
}

void inctc(char *cname,int npntsh,int *ncontrib,int **index,double **ctca)
{
 int i,j,idum;
 FILE *infile;
 double tmp;

 infile=fopen(cname,"r");
 fread(&idum,sizeof(int),1,infile); 
 if(idum!=npntsh) {printf("\nerror in file %s\n",cname); exit(1);} 
 for(i=1;i<=npntsh;i++)
    {
     fread(&ncontrib[i],sizeof(int),1,infile); 
     for(j=1;j<=ncontrib[i];j++) 
        {
         fread(&index[i][j],sizeof(int),1,infile);         	 
         fread(&tmp,sizeof(double),1,infile);         	 
         ctca[i][j]=tmp;
        }	 
    } 
 fclose(infile); 
}

double tauctcatau(int npntsh,int *ncontrib,int **index,double **ctca,double *tau)
{
 int i,j,k;
 double tmp,res;
 
 res=0.0;
 for(i=1;i<=npntsh;i++)
    {
     tmp=0;
     for(j=1;j<=ncontrib[i];j++)
        {
         k=index[i][j];
         tmp+=ctca[i][j]*tau[k];
        }
     res+=tau[i]*tmp;
    }
 return sqrt(res);
/* return res;*/
}


