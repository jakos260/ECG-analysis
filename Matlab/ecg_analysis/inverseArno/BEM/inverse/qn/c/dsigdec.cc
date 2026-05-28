/* file dsigdec.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

# define max(a,b) 		(a<b ? b : a)
# define min(a,b) 		(a>b ? b : a)

void dsigdec(int m,int n,double **u,double *sigma,double **v)
{
 int i,j;
 int iend,imin;
 double sigmin,tmp;

 iend=max(m,n);
/* if(iend>1) do*/
 if(iend>1) do
  {
   sigmin=sigma[iend];
   imin=iend;
   for (i=1;i<=iend;i++)
   {
    tmp=sigma[i];
    if (tmp<sigmin)
     {
      sigmin=tmp;
      imin=i;
     }
    }
   if (imin!=iend)
   {
    tmp=sigma[iend];          /*** exchange sigma-values ***/
    sigma[iend]=sigmin;
    sigma[imin]=tmp;
/*    printf("exchanged %d %d\n",iend,imin);*/
    for (i=1;i<=m;i++)
     {
      tmp=u[i][iend];
      u[i][iend]=u[i][imin];
      u[i][imin]=tmp;
      }
				/*** exchange u and v columns ***/
     for (i=1;i<=n;i++)
     {
      tmp=v[i][iend];
      v[i][iend]=v[i][imin];
      v[i][imin]=tmp;
      }
     }
     iend --;

    } while (iend>1);
   }


