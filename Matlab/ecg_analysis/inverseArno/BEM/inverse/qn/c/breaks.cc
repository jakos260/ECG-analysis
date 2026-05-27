/*****************************************************************************/
/* file breaks.c    Geertjan Huiskamp UC-Irvine June '92                     */
/*                 							     */
/* set up neighbor-definition int function neighbor.c                        */
/* determine breakthroughs    int function locsup.c                          */
/*									     */
/* external program should allocate *nnb **nb *supind *supval *supwidth      */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>

int neighbor(int npth,int ntrh,int **itrh,int *nnb,int **nb)
{
 int i,j,k,l;
 int ij,ik,maxnnb;

 printf("neighbor\n");
 for(i=1;i<=ntrh;i++) for(j=1;j<=3;j++)
    {
     k=j+1;
     if(k==4) k=1;
     ij=itrh[i][j];
     ik=itrh[i][k];
     if(ij>ik)
       {
        nnb[ij]+=1;
        nb[ij][nnb[ij]]=ik;
        nnb[ik]+=1;
        nb[ik][nnb[ik]]=ij;
       }
    }
 maxnnb=0; 
 for(i=1;i<=npth;i++) 
    {
     nb[i][0]=i; 
     if(nnb[i]>maxnnb) maxnnb++;
    }
 printf("maxnnb=%d\n",maxnnb);
 return maxnnb;
}

int locsup(double eps,int npth,int *nnb,int **nb,double *fun,
           int *supind,double *supval,double *supwidth)
{
 int   nsup,i,j,k,lomax,lomin;
 double cval,cvalp,cvalm,nbval,nbmean,nbmeans,width;

 nsup=0;
 for(i=1;i<=npth;i++)
    {
     cval=fun[i];
     cvalp=cval;
     cvalm=cval;
     lomin=1;
     lomax=1;
     nbmean=0;
     nbmeans=0;
     for(j=1;j<=nnb[i];j++)
        {
         nbval=fun[nb[i][j]];        
         if(nbval>cvalp) lomax=0; 
         if(nbval<cvalm) lomin=0;        
         nbmean+=nbval;
         nbmeans+=nbval*nbval;
        }       
     if(lomax || lomin)
       {
        width=cval-(nbmean/nnb[i]); 
        if(fabs(width)>=eps) 
          {
           ++nsup;
           supind[nsup]=i;
           supval[nsup]=cval;
           supwidth[nsup]=width;
          }
       }
    }
 return nsup;
}      
      

