/*  finddhk.c 

     vindt voor een electrode de dichtstbijzijnde driehoek
     door naar de maxdhk dichtstbijzijnden te kijken en genereert een
     electrode lambda-mu file met die driehoeken.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "trilib.h"
#define maxdhk 10000
#define herrie 0

double **A;
double *B;
double d;

static double det(double r1[3], double r2[3], double r3[3])
 {
  double det;

  det=r1[0]*(r2[1]*r3[2]-r2[2]*r3[1]);
  det=det-r1[1]*(r2[0]*r3[2]-r2[2]*r3[0]);
  det=det+r1[2]*(r2[0]*r3[1]-r2[1]*r3[0]);
  return det;
 }

extern void ludec(double **,int,double *);
extern void solvel(double **,int,double *);
extern void solveu(double **,int,double *);



double afstand2(double pnt1[3], double pnt2[3])
{
   double sum=0;
   int   i;
   for (i=0;i<3;i++)
     sum += (pnt1[i]-pnt2[i])*(pnt1[i]-pnt2[i]); 
   return(sqrt(sum));
}
   

void find_dhk2(double ppp[3],double c[3],double *la,double *mu,double *eta,
          int *eldhk, int npnt, double (*pnt)[3], int ndhk, int (*dhk)[3])

{
  int i,j,k,l,ready,mindhk[maxdhk],eldhk2=100,minmindhk;
  double s,sum,minsum[maxdhk];
  double a[3],b[3],p[3],d[3];
  double cond,deter;
  double mineta=1000,mineta2=1000,minla,minmu;
  double eta2=0,mu2=0,la2=0;
  int mmaxdhk=ndhk/2;

  A=matrix(1,3,1,3);
  B=vector(1,3);

  *la = 0;
  *mu = 1;
  *eldhk = 100;

  for (k=0;k<maxdhk;k++){
    minsum[k]=100;
    mindhk[k]=-1;
  }
  for (i=0;i<ndhk ;i++){
    sum =0;
    for (j=0;j<3;j++){
      sum += afstand2(ppp,pnt[dhk[i][j]]);
    }
    ready=0;
    for (j=0;(j<maxdhk)&&(!ready);j++){
      if (sum<minsum[j]){
        for (k=maxdhk-1;k>j;k--){
          minsum[k]=minsum[k-1];
          mindhk[k]=mindhk[k-1];
        }
        minsum[j]=sum;
        mindhk[j] = i;
        ready++;
      }
    } 
  }     
  ready = 0;
  for (k=0;(k<ndhk);k++){
    for (j=0;j<3;j++){  
        p[j]=pnt[dhk[mindhk[k]][0]][j];
        a[j]=pnt[dhk[mindhk[k]][1]][j]-pnt[dhk[mindhk[k]][0]][j];
        b[j]=pnt[dhk[mindhk[k]][2]][j]-pnt[dhk[mindhk[k]][0]][j];
    }
/*
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
*/
/*
    for (i=1;i<=3;i++){
      A[i][1] = a[i-1];
      A[i][2] = b[i-1];
      A[i][3] = c[i-1]; 
      B[i]    = ppp[i-1]-p[i-1];
    }
*/
    for (i=0;i<3;i++)
      d[i]=ppp[i]-p[i];
    deter=det(a,b,c);
    if (deter==0)
      continue;
    *la=det(d,b,c)/deter;
    *mu=det(a,d,c)/deter;
    *eta=det(a,b,d)/deter;
    if ((*la>=0)&&(*la<=1.0)&& (*mu>=0)&&(*mu+*la<=1.0)) {
      if (abs(*eta)<abs(mineta)){
        if (herrie)
          printf("dhk=%4d  la=%7.4f  mu=%7.4f  eta=%9.4f\n",
                                      mindhk[k],*la,*mu,-*eta);
        mineta = -*eta;
        minmindhk = mindhk[k];
        minla = *la;
        minmu = *mu;
        ready=1;
      }
    }
  }
  *la=minla;
  *mu=minmu;
  *eta=mineta;
  *eldhk=minmindhk;          
 }
/*
    ludec(A,3,&cond);
    solvel(A,3,B);
    solveu(A,3,B);
    if ((B[1]>=0)&&(B[1]<=1.0)&& (B[2]>=0)&&(B[2]<=1.0)&&((B[1]+B[2])<=1.0)){
      if (abs(B[3])<abs(mineta)){
        if (herrie)
          printf("dhk=%4d  la=%7.4f  mu=%7.4f  eta=%9.4f\n",
                                      mindhk[k],B[1],B[2],-B[3]);
        mineta = -B[3];
        *eldhk = mindhk[k];
        *la = B[1];
        *mu = B[2];
        ready=1;
        *eta = -B[3];
      }
    }
    else{
      if (abs(B[3])<abs(mineta2)){
        mineta2=-B[3];
        eldhk2=mindhk[k];
        la2=B[1];
        mu2=B[2];
        eta2=-B[3];
      }
    }
  } 
*/
/*printf("%7.4f %7.4f\n",*eta,eta2);*/
/*
  if ((!ready)){
    la2=(la2<0)?0:(la2>1)?1:la2;
    mu2=(mu2<0)?0:(mu2>1)?1:mu2;
    if ((s=la2+mu2) > 1){
     la2 /=s;
     mu2 /=s; 
    }
    if (herrie){
      printf("Didn't succeed in finding projection. Using best I have:\n");
      printf("dhk=%4d  la=%7.4f  mu=%7.4f  eta=%9.4f\n", eldhk2,la2,mu2,eta2);
    }
    *eldhk =eldhk2;
    minla = la2;
    minmu = mu2;
    mineta = eta2;
  }
}
*/

