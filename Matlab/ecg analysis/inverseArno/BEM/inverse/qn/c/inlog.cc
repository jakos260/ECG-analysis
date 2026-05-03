/* inlog.c initial estimate lcurve generation */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include "trilib.h"

#define EPS 1.e-07
#define BIG 1e+34

#define max(a,b) 		(a<b ? b : a)
#define min(a,b) 		(a>b ? b : a)

/*
extern double *vector(int,int);
extern int   *ivector(int,int);
extern double **getmatrix(char *,int *,int *);
extern double **matrix(int,int,int,int);
extern int   **imatrix(int,int,int,int);
extern void  free_matrix(double **,int,int,int,int);
*/

extern int   neighbor(int,int,int **,int *,int **);
extern void  svdgj(double **,int,int,double *,double **);
extern void  sigdec(int,int,double **,double *,double **);
extern double **vcalcu(int,int,int,int,double,double *,double **);
extern int   locsup(double,int,int *,int **,double *,int *,double *,double *);

double *rimp(int,int,double **);
void  inctc(char *,int ,int *,int **,double **);
void  trpsinv(int,int,int,double **,double *,double **,double **);
void  calctau(int,int,double **,double *,double *,double *,double *);
int   shift(int,int,double *,double *,double *,int *);
double tauctcatau(int,int *,int **,double **,double *);
int   optims(int,int,int,double **,double **,double,double *);
double referr(int,int,int,double *,double *);

int main(int argc, char *argv[])
{
 int i,j,npntsh,ntrsh,idum,maxnnb,nel,iwind,itref,iwidth,irank,ibeg,iend;
 int nat,ishift,nsup,noref,lshift,nnbmax=50;
 int **itrh,*nnb,**nb,*ncontrib,**cindex,*supind;
 double **a,**u,**v,*sig,**vm,*vin,**ctca,*tau;
 double **psi,**vi,*tref,*supval,*supwidth;
 double fdum,tmp,src,str,resn,treg,rest,tres,sigmax,tmin,tmax,norm;
 char heartname[80],amaname[80],bsmname[80],refname[80],ctcname[80];
 char line[80];
 FILE *fp,*flog;

 setbuf(stdout,NULL);

 printf("triangulation heart>");
 gets(heartname);
 fp=fopen(heartname,"r");
 fscanf(fp,"%d\n",&npntsh);
 for(i=1;i<=npntsh;i++) 
 fscanf(fp,"%d %lf %lf %lf\n",&idum,&fdum,&fdum,&fdum);
 fscanf(fp,"%d\n",&ntrsh);
 itrh=imatrix(1,ntrsh,1,3);
 for(i=1;i<=ntrsh;i++) 
 fscanf(fp,"%d %d %d %d\n",&idum,&itrh[i][1],&itrh[i][2],&itrh[i][3]);
 fclose(fp);
 printf("file read\n");

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
 cindex=imatrix(1,npntsh,1,nnbmax);
 ncontrib=ivector(1,npntsh);
 ctca=matrix(1,npntsh,1,nnbmax);
 inctc(ctcname,npntsh,ncontrib,cindex,ctca);

 u=matrix(1,nel,1,npntsh); 
 v=matrix(1,npntsh,1,npntsh);
 sig=vector(1,npntsh);
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
 psi=matrix(1,npntsh,1,nel);
 
 ibeg=1;
 iend=min(nel,npntsh);

 irank=ibeg;
 
 iwidth=1;
 flog=fopen("lcurve.tab","w");
 fprintf(flog,"%d 8\n",iend-1);
 printf ("\nran  smax  tmin  tmax width      reg     res    tres    desc nsup\n");

 while (irank<iend && sig[irank]>=EPS)	 
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
   lshift=shift(npntsh,itref,tau,&tmin,&tmax,&iwidth);
     {
      if(lshift)
        {
         treg=tauctcatau(npntsh,ncontrib,cindex,ctca,tau);
         vi=vcalcu(nel,npntsh,iwidth,iwind,str,tau,a);    
         ishift=optims(nel,itref,iwidth,vm,vi,resn,&rest);
         free_matrix(vi,1,nel,1,iwidth);  
         tres=referr(noref,npntsh,ishift,tau,tref);     
         nsup=locsup(1.0,npntsh,nnb,nb,tau,supind,supval,supwidth);      
        }
      if(ibeg!=iend)
        {
          fprintf(flog,"%3d\t%6.4f\t%6.4f\t%8.2f\t%6.4f\t%8.4f\t%8.4f\t%d\n",
	               irank,sig[irank]/sigmax,norm,treg,rest,tres,fdum,nsup);
          printf("\n%3d %5.3f %5.2f %5.2f   %3d%10.2f%8.4f%8.4f%8.4f  %2d",
          irank,sig[irank]/sigmax,tmin+ishift,tmax+ishift,iwidth,treg,rest,tres,fdum,nsup);
        }     
      }
   irank++;
  }               
 printf ("\n");
 fclose (flog);
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

int shift(int npntsh,int itref,double *tau,double *tmin,double *tmax,int *iwidth)
{
 int i;
 double width;

 for(i=1;i<=npntsh;i++) tau[i]-=(*tmin)-1;
 *tmax-=(*tmin)-1;
 *tmin=1;
 width=*tmax-*tmin;
 *iwidth=((int)(width))+1;
 if(*iwidth<itref)
   return 1;
 else
   {
    printf("*");
    return 0;
   }
}

int optims(int nel,int itref,int iwidth,
           double **vm,double **vi,double resn,double *rest)
{
 int is,ie,it,iw,ishift;
 double resmin,res,vv,ww;

 resmin=BIG;
 for (is=1;is<=(itref-iwidth);is++)
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

void inctc(char *cname,int npntsh,int *ncontrib,int **cindex,double **ctca)
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
         fread(&cindex[i][j],sizeof(int),1,infile);         	 
         fread(&tmp,sizeof(double),1,infile);         	 
         ctca[i][j]=tmp;
        }	 
    } 
 fclose(infile); 
}

double tauctcatau(int npntsh,int *ncontrib,int **cindex,double **ctca,double *tau)
{
 int i,j,k;
 double tmp,res;
 
 res=0.0;
 for(i=1;i<=npntsh;i++)
    {
     tmp=0;
     for(j=1;j<=ncontrib[i];j++)
        {
         k=cindex[i][j];
         tmp+=ctca[i][j]*tau[k];
        }
     res+=tau[i]*tmp;
    }
 return sqrt(res);
/* return res;*/
}


