/***************************************************************************/
/* File allfuns.c                                                           */
/* Maria Orgonasova Geertjan Huiskamp                                      */
/* Medical physics and biophysics                                          */
/* University of Nijmegen The Netherlands, November 1988 (.c december l991)*/
/*                                                                         */
/* Minimalization-function definition for program qnall                    */
/* Evaluation of function fun=funall(), and grad[gradient in fugrad () with     */
/* respect to parameters tau[j] all leads, all times.                      */
/* Simulation preformed with quadratic activation function actfun() with   */
/* window of extend iwind.                                                 */
/*                                                                         */
/* Function qnscl ()                                                       */
/* Calculation of balance between error and regularisation in order to be  */
/* able to use a 'normalized' regularization parameter basically a         */
/* simplified funall() function                                            */
/*                                                                         */
/* function funmon() monitoring function                                   */
/***************************************************************************/
    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*#include "matg.h"*/

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define MIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) < (maxarg2) ?\
	(maxarg1) : (maxarg2))

extern double *grad;
extern double **at,**vs;

extern double **ctca;

extern double **a,**v;
/*extern double **a,**v,**ctca;*/
extern int   *ncontrib,**cindex,*itmi,*itpl;
extern int   nel,iwind,jwind,npntsh,nat,itref,locext;
extern double regmin,regmax,regstep,regpar,err,vrms,trefn,tmin,tmax,upper,lower;
extern double reg;
extern double *row,*tref,*bounds,terror,epsbreak;
extern char  qref[20],bon[1],lis[1];
extern FILE  *tlisten;

int          ie,itm,itp,it,nfc=0;
double       ayx,taux;

double funall(double *);
void  fugrad(double *,double *);
int   funmon(FILE *,int,double *,double *);
double qnscl(double *);
void  conrp(double *);
double actfun(int,double,double);
void  actval(double *);
void  crevsd(double *);
void  integvs(void);
void  reserr(void);

extern double *supval,*supwidth;
extern int *nnb,**nb,*supind;
int locsup(double,int,int *,int **,double *,int *,double *,double *);

double funall (double *tau)
/**************************/
{
 int i;
 double fun;

 ++nfc;
 err=0;
 reg=0;
 conrp(tau);
 actval(tau);
 for(ie=1;ie<=nel;ie++)
    {
     crevsd(tau);
     integvs(); 
     reserr();
    }
 fun=err+reg*regpar;
 return (fun);
}


double qnscl(double *tau)
{
 double bal;

 funall(tau);
 nfc=0;
 printf("err=%g , reg=%g \n",err,sqrt(fabs(reg)));
 bal=err/reg;
 printf("balance=%e\n",bal);
 return bal;
}

void fugrad(double *tau,double *fgrad)     
/***************** !! must be called after funall !! **/
{
 int i,iel,k,itm1;
 double tmp;

 for(i=1;i<=npntsh;i++) fgrad[i]=grad[i];
 for(iel=1;iel<=nel;iel++)
    for(i=1;i<=npntsh;i++)
       {
        tmp=0;
        ayx=2*a[iel][i];
        itm=itmi[i];
        itm1=itm-1;
        itp=itpl[i];
        for(k=itm;k<=itp;k++)
           {
            it=k-itm1;
            tmp+=vs[iel][k]*at[i][it]*ayx;
           }
        fgrad[i]-=tmp;
       }
}

void conrp(double *tau)
/*********************/
{
 int i,j,k;
 double tmp;
 
 for(i=1;i<=npntsh;i++)
    {
     tmp=0;
     for(j=1;j<=ncontrib[i];j++)
        {
         k=cindex[i][j];
         tmp+=ctca[i][j]*tau[k];
        }
     reg+=tau[i]*tmp;
     grad[i]=2*regpar*tmp;
    }
}

void actval(double *tau)
/***********************/
{
 int i,k,itm1;
 double t;

 for(i=1;i<=npntsh;i++)
    {
     taux=MAX(MIN(tau[i],upper),lower);
     if(bon[0]=='y') bounds[i]=taux;
     itmi[i]=itm=taux-jwind+1;
     itm1=itm-1;
     itpl[i]=itp=taux+jwind;
     for(k=itm;k<=itp;k++)
        {
         t=(double)k;
         it=k-itm1;
         at[i][it]=actfun(iwind,t,taux);
        }
    }
}

void crevsd(double *tau)
/***********************/
{
 int i,j,k,itm1;

 for(j=1;j<=itref;j++) row[j]=0;
 for(i=1;i<=npntsh;i++)
    {
     ayx=a[ie][i];
     itm=itmi[i];
     itm1=itm-1;
     itp=itpl[i];
     for(k=itm;k<=itp;k++)
        {
         it=k-itm1;
         row[k]+=at[i][it]*ayx;
        }
    }
}

void integvs(void)
/*******************/
{
 int j,j1;

 vs[ie][1]=row[1]/2;
 for (j=2;j<=itref;j++) { j1=j-1; vs[ie][j]=vs[ie][j1]+(row[j]+row[j1])/2;}
}

void reserr(void)
/******************/
{
 int j;
 double diff;
 FILE *vv;

 for(j=1;j<=itref;j++)
    {
     diff=vs[ie][j]-v[ie][j];
     err+=diff*diff;
     vs[ie][j]=diff;
    }
}

double actfun(int iwind,double t,double tau)
/********************************************/
{
 double hw,dt;

 hw=iwind/2;
 dt=t-tau;
 if(dt<0) return ((1+dt/hw)/hw); else return ((1-dt/hw)/hw);
}

int funmon(FILE *mon,int nit,double *p,double *g)
/**************************************************/
{
 double res,regn,tres,refi,taui,ngrad;
 int i,nsup; 

 res=sqrt(err)/vrms;    
 regn=sqrt(reg);
 tres=0;
 ngrad=0;
 if(qref[0]=='*') trefn=0;
 if(bon[0]=='y') for(i=1;i<=npntsh;i++) p[i]=bounds[i];
 if(lis[0]=='y') 
   {
    tlisten=fopen("listen.tim","w"); 
    fprintf(tlisten,"%4d\n",npntsh);
    for(i=1;i<=npntsh;i++) fprintf(tlisten,"%.2f\n",p[i]);   
    fclose(tlisten);
   }
 tmin=nat;
 tmax=0;
 for(i=1;i<=npntsh;i++)
    {
     refi=tref[i];
     taui=p[i];      
     tmax=MAX(taui,tmax);
     tmin=MIN(taui,tmin);
     if(qref[0]=='*')
       {
        tref[i]=taui;
        trefn+=tref[i]*tref[i];
       }
     refi-=taui;
     tres+=refi*refi;
     ngrad+=g[i]*g[i];
    }
 terror=tres=sqrt(tres/trefn); 
 ngrad=sqrt(ngrad); 
 locext=nsup=locsup(epsbreak,npntsh,nnb,nb,p,supind,supval,supwidth);
 res=sqrt(err)/vrms;    
 fprintf(mon,"%d\t %d\t %.4f\t %6.0f\t %.4f\t %e\t %d\n",
              nit,nfc,res,regn,ngrad,tres,nsup); 
 if(lis[0]=='y') 
   printf("%4d %4d %3.1f\t %3.1f\t %.4f\t %6.0f\t %.4f\t %e \t %d\n",
             nit,nfc,tmin,tmax,res,regn,ngrad,tres,nsup); 
 nfc=0;
 if (tres==0)
   return true;
 else
   return false;
}
