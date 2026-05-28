
  /**********************************************************************/
  /*                                                                    */
  /*    parmelec.c                                                      */
  /*                                                                    */
  /*    Program to determine the "lambda-mu" parameters	of the		*/
  /*	BSM-electrodes based on the PARMA system on a triangulated torso. */
  /*									*/
  /*	The assignment of electrodes closely resembles the way it is	*/
  /*	done by BSM laborants.						*/
  /*									*/
  /*    Rudi Hoekema, Feb 22, 1195, after program eegelec.c by		*/
  /*									*/
  /*    Thom Oostendorp							*/
  /*    November 25, 1993						*/
  /*									*/
  /**********************************************************************/
			    
/*
        input:
	      - standard triangle file defining the geometry of the torso 

        output:
	      - electrode description file
                (file containing the "lambda-mu" parameters of the electrodes).
										
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "trilib.h"

#define NEL	    120
#define distance     0.05

int 	ndhk, npnt;
int 	(*dhk)[3]=NULL;
double 	(*pnt)[3]=NULL, dir[3],p[3],p1[3],p2[3],p3[3],n[3];
double   pl1[3],pr1[3],pf1[3],pl2[3],pr2[3],pf2[3];
double	(*refpnt)[3];
int	nel;
CONTOUR tmp[NEL], tmp2[NEL], tmp1[NEL], elec[NEL],el[NEL];
double   norm_x[3]={1,0,0};
double   norm_minx[3]={-1,0,0};
double   norm_y[3]={0,1,0};
double   norm_miny[3]={0,-1,0};
double   norm_z[3]={0,0,1};
double   norm_xy1[3]={1,1,0};
double   norm_xy2[3]={-1,1,0};
double   norm_minz[3]={0,0,-1};
double   origin[3]={0,0,0};
double   vertshift=0.0;



char elname[NEL][21]=
        {
           "1", "2","3","4","5","6","7","8","9",
           "10","11","12","13","14","15","16","17","18","19",
           "20","21","22","23","24","25","26","27","28","29",
           "30","31","32","33","34","35","36","37","38","39",
           "40","41","42","43","44","45","46","47","48","49",
           "50","51","52","53","54","55","56","57","58","59",
           "60","61","62","63","64","65","66","67","68","69",
           "70","71","72","73","74","75","76","77","78","79",
           "80","81","82","83","84","85","86","87","88","89",
           "90","91","92","93","94","95","96","97","98","99",
           "100","101","102","103","104","105","106","107","108","109",
           "110","111","112","113","114","115","116","117","118","119",
           "120"
        };

double afstand(double pnt1[3], double pnt2[3]);
void find_dhk(double ppp[3],double *la,double *mu,double *eta,int *eldhk,
        int npnt, double (*pnt)[3], int ndhk, int (*dhk)[3]) ;

int assign_el(char *name, CONTOUR *c, int i, CONTOUR *el)
 {
  int j;
  
  for (j=0; j<NEL; j++)
    if (strcmp(elname[j], name)==0)
      break;
  if (j==NEL)
   {
    printf("Electrode %s undefined\n", name);
    return -1;
   }
  el[j].tr=c[i].tr;
  el[j].l=c[i].l;
  el[j].m=c[i].m;
  return j;
 }

int assign_lm_el(char *name, int tr, double l, double m, CONTOUR *el)
 {
  int j;
  
  for (j=0; j<NEL; j++)
    if (strcmp(elname[j], name)==0)
      break;
  if (j==NEL)
   {
    printf("Electrode %s undefined\n", name);
    return -1;
   }
  el[j].tr=tr;
  el[j].l=l;
  el[j].m=m;
  return j;
 }

double cont_len(int nlm,DCONT *c,
            int npnt, double (*pnt)[3], int ndhk, int (*dhk)[3])
{
  double som=0,r1[3],r2[3];
  int i;

  for (i=0;i<nlm;i++){
    routlm(r1, c[i].tr, c[i].l, c[i].m, npnt, pnt, ndhk, dhk);
    routlm(r2, c[i].tr, c[i].l2, c[i].m2, npnt, pnt, ndhk, dhk);
    som += afstand(r1,r2);
  }
  return som;
}

  

int compose_normal(double n[3],int tr1,int tr2,
                  int npnt, double (*pnt)[3], int ndhk, int (*dhk)[3])
{
  double n1[3],n2[3],a[3],b[3],len1,len2;
  int i;

  for (i=0;i<3;i++){
    p[i]=pnt[dhk[tr1][0]][i];
    a[i]=pnt[dhk[tr1][1]][i]-pnt[dhk[tr1][0]][i];
    b[i]=pnt[dhk[tr1][2]][i]-pnt[dhk[tr1][0]][i];
  }
  len1 = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  len2 = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  for (i=0;i<3;i++){
    a[i]/=len1;
    b[i]/=len2;
  }
  n1[0]=a[1]*b[2]-a[2]*b[1];
  n1[1]=a[2]*b[0]-a[0]*b[2];
  n1[2]=a[0]*b[1]-a[1]*b[0];

  for (i=0;i<3;i++){
    p[i]=pnt[dhk[tr2][0]][i];
    a[i]=pnt[dhk[tr2][1]][i]-pnt[dhk[tr2][0]][i];
    b[i]=pnt[dhk[tr2][2]][i]-pnt[dhk[tr2][0]][i];
  }
  len1 = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  len2 = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  for (i=0;i<3;i++){
    a[i]/=len1;
    b[i]/=len2;
  }
  n2[0]=a[1]*b[2]-a[2]*b[1];
  n2[1]=a[2]*b[0]-a[0]*b[2];
  n2[2]=a[0]*b[1]-a[1]*b[0];

  for (i=0;i<3;i++)
    n[i]=0.5*(n1[i]+n2[i]);
}



int contcross(DCONT *c1, int f1, int l1, DCONT *c2, int f2, int l2,
		int *tr, double *l, double *m, double criterium,
		int npnt, double (*pnt)[3], int ndhk, int (*dhk)[3])
 {
  int i, j, k, tra, trb;
  double p1[3], p2[3], s1[3], s2[3], min=1e20;
  double la, ma, lb, mb, r1, r2;
  
  for (i=f1; i<l1; i++)
   {
    routlm(p1, c1[i].tr, c1[i].l, c1[i].m, npnt, pnt, ndhk, dhk);
    routlm(s1, c1[i].tr, c1[i].l2, c1[i].m2, npnt, pnt, ndhk, dhk);
    for (k=0; k<3; k++)
      s1[k] -= p1[k];
    for (j=f2; j<l2; j++)
     {
      routlm(p2, c2[j].tr, c2[j].l, c2[j].m, npnt, pnt, ndhk, dhk);
      routlm(s2, c2[j].tr, c2[j].l2, c2[j].m2, npnt, pnt, ndhk, dhk);
      for (k=0; k<3; k++)
        s2[k] -= p2[k];
      if (linesect(p1, s1, p2, s2, &r1, &r2))
       {
        if (r1>=0 && r1<=1 && r2>=0 && r2<=1)
	 {
	  tra=c1[i].tr;
	  la=(1-r1)*c1[i].l+r1*c1[i].l2;
	  ma=(1-r1)*c1[i].m+r1*c1[i].m2;
	  trb=c2[j].tr;
	  lb=(1-r2)*c2[j].l+r2*c2[j].l2;
	  mb=(1-r2)*c2[j].m+r2*c2[j].m2;
	  routlm(p2, tra, la, ma, npnt, pnt, ndhk, dhk);
	  routlm(s2, trb, lb, mb, npnt, pnt, ndhk, dhk);
	  if (vecdist(p2, s2)<=criterium)
	   {
	    *tr=tra;
	    *l=la;
	    *m=ma;
	    return 1;
	   }
	 }
       }
     }
   }
  return 0;
 }

void write_elec(char *filename)
 {
  int	i;
  FILE	*file;

  nel=NEL;
  file=fopen(filename, "w");
  if (file==NULL)
   {
    printf("Error opening %s\n", filename);
    exit(1);
   }
  fprintf(file, "%4d\n", nel);
  for (i=0; i<nel; i++)
    if (fprintf(file,"%4d %7.4f %7.4f !%s\n", elec[i].tr+1, 
		elec[i].l, elec[i].m, elname[i])==EOF)
     {
      printf("Error writing to %s\n", filename);
      exit(1);
     }
  fclose(file);
 }


void out_of_mem(void)
 {
  printf("Sorry, out of memory");
  exit(1);
 }


void main(int argc, char *argv[])
 {
  int	i, j, k, tr,pp,ok;
  int	nlmh1,nlmh2, index;
  DCONT *h1,*h2;
  double l, m;
  char 	filename[80]="\0";
  FILE	*file;
  double pos[38],pos2[38];
  double ymin,ymax;
  int cmin,cmax;
  DCONT *c,*vc;
  int nlm,vnlm;
  char line[80]="\0";


/* get user input */

  if (argc>1)
    strcpy(filename, argv[1]);
  else
   {
    printf("Name triangulated torso: ");
    gets(filename);
   }
  if (filename[0]=='\0')
    exit(0);
  if (!dhkinp(filename, &npnt, &ndhk, &pnt, &dhk))
    exit(1);

  if (argc>2)
    strcpy(filename, argv[2]);
  else
   {
    printf("Name output file: ");
    gets(filename);
   }
  if (filename[0]=='\0')
    exit(0);
  
  if (argc>3)
    strcpy(line, argv[3]);
  else if (argc==1)
   {
    printf("Shift points vertically: [m] ");
    gets(line);
   }
  if (strlen(line)){
    sscanf(line,"%lf",&vertshift);
  }


  h1=calloc(ndhk, sizeof(DCONT));
  if (h1==NULL)
    out_of_mem();

  h2=calloc(ndhk, sizeof(DCONT));
  if (h2==NULL)
    out_of_mem();

  c=calloc(ndhk, sizeof(DCONT));
  if (c==NULL)
    out_of_mem();

  vc=calloc(ndhk, sizeof(DCONT));
  if (vc==NULL)
    out_of_mem();

  for (i=0;i<18;i++)
    pos[i]=i/9.0;
  
  for (i=0;i<18;i++)
    pos2[i]=i*distance;

/*  get wilson's electrode positions */

  p1[0]=0; p1[1]= 1; p1[2]=0.125+vertshift;
  contour(p1, norm_y, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  routlm(pl1,h1[0].tr,h1[0].l,h1[0].m,npnt,pnt,ndhk,dhk);

  p1[0]=0; p1[1]=-1; p1[2]=0.125+vertshift;
  contour(p1, norm_y, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  routlm(pr1,h1[0].tr,h1[0].l,h1[0].m,npnt,pnt,ndhk,dhk);
  
  contour(pr1, pl1, 1, norm_x,0, norm_y, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);


  elecs(nlmh1, h1, 0, index, npnt, pnt, ndhk, dhk, 10, pos, 1, tmp);
  assign_el("1", tmp, 0, elec);
  assign_el("2", tmp, 9, elec);

  p1[0]=0; p1[1]=1; p1[2]=-0.25+vertshift;
  ok=0; 
  for (i=0;(i<npnt) && !ok;i++){
    if ((pnt[i][0]==0)&&(pnt[i][2]==p1[2])&&(pnt[i][1]>0)){
      for (j=0;j<3;j++)
        pl1[j]=pnt[i][j];
      ok++;
    }
  }
  if (!ok){
    contour(p1, norm_y, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
      &nlmh2, h2, &index);
    routlm(pl1,h2[0].tr,h2[0].l,h2[0].m,npnt,pnt,ndhk,dhk);
  }

  contour(pl1, norm_y,0,  norm_x,0, norm_y, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  elecs(nlmh1, h1, 0, index, npnt, pnt, ndhk, dhk, 10, pos, 1, tmp);
  assign_el("3", tmp, 0, elec);

/*  get reference points at sternum */

  p1[0]= 1.0; p1[1]= 0.0; p1[2]= 0.0+vertshift;
  ok=0; 
  for (i=0;(i<npnt) && !ok;i++){
    if ((pnt[i][1]==0)&&(pnt[i][2]==p1[2])&&(pnt[i][0]>0)){
      for (j=0;j<3;j++)
        pf1[j]=pnt[i][j];
      ok++;
    }
  }
  if (!ok){
    contour(p1, norm_x, 0, norm_y,0, norm_y, npnt, pnt, ndhk, dhk, 
      &nlmh1, h1, &index);
    routlm(pf1,h1[0].tr,h1[0].l,h1[0].m,npnt,pnt,ndhk,dhk);
  }

  p1[0]= 1.0; p1[1]= 0.0; p1[2]= -0.15+vertshift;
  contour(p1, norm_x, 0, norm_y,0, norm_y, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  routlm(pf2,h1[0].tr,h1[0].l,h1[0].m,npnt,pnt,ndhk,dhk);

/*  get reference points at left side */

  p1[0]=0; p1[1]=1; p1[2]=0+vertshift;
 
  ok=0; 
  for (i=0;(i<npnt) && !ok;i++){
    if ((pnt[i][0]==0)&&(pnt[i][2]==p1[2])&&(pnt[i][1]>0)){
      for (j=0;j<3;j++)
        pl1[j]=pnt[i][j];
      ok++;
    }
  }
  if (!ok){
    contour(p1, norm_y, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
      &nlmh2, h2, &index);
    routlm(pl1,h2[0].tr,h2[0].l,h2[0].m,npnt,pnt,ndhk,dhk);
  }

  p1[0]=0; p1[1]=1; p1[2]=-0.15+vertshift;
  contour(p1, norm_y, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
    &nlmh2, h2, &index);
  routlm(pl2,h2[0].tr,h2[0].l,h2[0].m,npnt,pnt,ndhk,dhk);


/*  get reference points at right side */

  p1[0]=0; p1[1]=-1; p1[2]=0+vertshift;
  ok=0; 
  for (i=0;(i<npnt) && !ok;i++){
    if ((pnt[i][0]==0)&&(pnt[i][2]==p1[2])&&(pnt[i][1]<0)){
      for (j=0;j<3;j++)
        pr1[j]=pnt[i][j];
      ok++;
    }
  }
  if (!ok){
    contour(p1, norm_miny, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
      &nlmh2, h2, &index);
    routlm(pr1,h2[0].tr,h2[0].l,h2[0].m,npnt,pnt,ndhk,dhk);
  }

  p1[0]=0; p1[1]=-1; p1[2]=-0.15+vertshift;
  contour(p1, norm_miny, 0, norm_x,0, norm_x, npnt, pnt, ndhk, dhk, 
    &nlmh2, h2, &index);
  routlm(pr2,h2[0].tr,h2[0].l,h2[0].m,npnt,pnt,ndhk,dhk);

/*
printf("%7.4f %7.4f %7.4f\n",pf1[0],pf1[1],pf1[2]);
printf("%7.4f %7.4f %7.4f\n",pf2[0],pf2[1],pf2[2]);
printf("%7.4f %7.4f %7.4f\n",pl1[0],pl1[1],pl1[2]);
printf("%7.4f %7.4f %7.4f\n",pl2[0],pl2[1],pl2[2]);
printf("%7.4f %7.4f %7.4f\n",pr1[0],pr1[1],pr1[2]);
printf("%7.4f %7.4f %7.4f\n",pr2[0],pr2[1],pr2[2]);
*/
  /* strips on back side : construct horizontals */


  for (i=0;i<38;i++)
    pos[i]=i/32.0;
  
  contour(pl1, pr1, 1, norm_minx,0, norm_minx, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  elecs(nlmh1, h1, 0, index, npnt, pnt, ndhk, dhk, 33, pos, 1, tmp1);
  contour(pl2, pr2, 1, norm_minx,0, norm_minx, npnt, pnt, ndhk, dhk, 
    &nlmh2, h2, &index);
  elecs(nlmh2, h2, 0, index, npnt, pnt, ndhk, dhk, 33, pos, 1, tmp2);
  

  /* strip K */

  pp = 0;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("76",el,0,elec);
  assign_el("77",el,1,elec);
  assign_el("78",el,2,elec);
  assign_el("79",el,3,elec);



  /* strip L */

  pp = 3;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("80",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("81",el,0,elec);
  assign_el("82",el,1,elec);
  assign_el("83",el,2,elec);
  assign_el("84",el,3,elec);

  /* strip M */

  pp = 7;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("85",el,2,elec);
  assign_el("86",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("87",el,0,elec);
  assign_el("88",el,1,elec);
  assign_el("89",el,2,elec);
  assign_el("90",el,3,elec);

  /* strip N */

  pp = 12;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("91",el,4,elec);
  assign_el("92",el,3,elec);
  assign_el("93",el,2,elec);
  assign_el("94",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("95",el,0,elec);
  assign_el("96",el,1,elec);
  assign_el("97",el,2,elec);
  assign_el("98",el,3,elec);

  /* strip O */

  pp = 16;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("99",el,4,elec);
  assign_el("100",el,3,elec);
  assign_el("101",el,2,elec);
  assign_el("102",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("103",el,0,elec);
  assign_el("104",el,1,elec);
  assign_el("105",el,2,elec);
  assign_el("106",el,3,elec);

  /* strip P */

  pp = 20;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("107",el,4,elec);
  assign_el("108",el,3,elec);
  assign_el("109",el,2,elec);
  assign_el("110",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("111",el,0,elec);
  assign_el("112",el,1,elec);
  assign_el("113",el,2,elec);
  assign_el("114",el,3,elec);

  /* strip P */

  pp = 25;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("115",el,2,elec);
  assign_el("116",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("117",el,0,elec);
  assign_el("118",el,1,elec);
  assign_el("119",el,2,elec);
  assign_el("120",el,3,elec);


  /* strip P */

  pp = 32;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("4",el,0,elec);
  assign_el("5",el,1,elec);
  assign_el("6",el,2,elec);
  assign_el("7",el,3,elec);


  /* construct right horizontal */

  for (i=0;i<18;i++)
    pos[i]=i/16.0;

  contour(pf1, pr1, 1, norm_minx,0, norm_miny, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  elecs(nlmh1, h1, 0, index, npnt, pnt, ndhk, dhk, 17, pos, 1, tmp1);

  contour(pf2, pr2, 1, norm_minx,0, norm_miny, npnt, pnt, ndhk, dhk, 
    &nlmh2, h2, &index);
  elecs(nlmh2, h2, 0, index, npnt, pnt, ndhk, dhk, 17, pos, 1, tmp2);


  /* strip B */

  pp = 10;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("8",el,3,elec);
  assign_el("9",el,2,elec);
  assign_el("10",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("11",el,0,elec);
  assign_el("12",el,1,elec);
  assign_el("13",el,2,elec);
  assign_el("14",el,3,elec);

  /* strip C */

  pp = 6;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("15",el,3,elec);
  assign_el("16",el,2,elec);
  assign_el("17",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("18",el,0,elec);
  assign_el("19",el,1,elec);
  assign_el("20",el,2,elec);
  assign_el("21",el,3,elec);

  /* strip D */

  pp = 2;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("22",el,3,elec);
  assign_el("23",el,2,elec);
  assign_el("24",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("25",el,0,elec);
  assign_el("26",el,1,elec);
  assign_el("27",el,2,elec);
  assign_el("28",el,3,elec);

  /* strip E */

  pp = 0;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("29",el,4,elec);
  assign_el("30",el,3,elec);
  assign_el("31",el,2,elec);
  assign_el("32",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("33",el,0,elec);
  assign_el("34",el,1,elec);
  assign_el("35",el,2,elec);
  assign_el("36",el,3,elec);

  /* construct left horizontal */

  contour(pf1, pl1, 1, norm_minx,0, norm_y, npnt, pnt, ndhk, dhk, 
    &nlmh1, h1, &index);
  elecs(nlmh1, h1, 0, index, npnt, pnt, ndhk, dhk, 17, pos, 1, tmp1);

  contour(pf2, pl2, 1, norm_minx,0, norm_y, npnt, pnt, ndhk, dhk, 
    &nlmh2, h2, &index);
  elecs(nlmh2, h2, 0, index, npnt, pnt, ndhk, dhk, 17, pos, 1, tmp2);


  /* strip F */

  pp = 2;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("37",el,3,elec);
  assign_el("38",el,2,elec);
  assign_el("39",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("40",el,0,elec);
  assign_el("41",el,1,elec);
  assign_el("42",el,2,elec);
  assign_el("43",el,3,elec);

  /* strip G */

  pp = 4;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("44",el,3,elec);
  assign_el("45",el,2,elec);
  assign_el("46",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("47",el,0,elec);
  assign_el("48",el,1,elec);
  assign_el("49",el,2,elec);
  assign_el("50",el,3,elec);

  /* strip H */

  pp = 6;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("51",el,3,elec);
  assign_el("52",el,2,elec);
  assign_el("53",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("54",el,0,elec);
  assign_el("55",el,1,elec);
  assign_el("56",el,2,elec);
  assign_el("57",el,3,elec);

  /* strip I */

  pp = 8;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("58",el,3,elec);
  assign_el("59",el,2,elec);
  assign_el("60",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("61",el,0,elec);
  assign_el("62",el,1,elec);
  assign_el("63",el,2,elec);
  assign_el("64",el,3,elec);

  /* strip J */

  pp = 10;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_z, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("65",el,3,elec);
  assign_el("66",el,2,elec);
  assign_el("67",el,1,elec);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("68",el,0,elec);
  assign_el("69",el,1,elec);
  assign_el("70",el,2,elec);
  assign_el("71",el,3,elec);

  /* strip K+ */

  pp = 13;
  routlm(p1,tmp1[pp].tr,tmp1[pp].l,tmp1[pp].m, npnt, pnt, ndhk, dhk);
  routlm(p2,tmp2[pp].tr,tmp2[pp].l,tmp2[pp].m, npnt, pnt, ndhk, dhk);
  compose_normal(n,tmp1[pp].tr, tmp2[pp].tr,npnt,pnt,ndhk,dhk);

  contour(p1, n,0,p2, 1, norm_minz, npnt, pnt, ndhk, dhk, &vnlm, vc, &index);
  elecs(vnlm, vc, 0, 7, npnt, pnt, ndhk, dhk, 17, pos2, 0, el);
  assign_el("72",el,0,elec);
  assign_el("73",el,1,elec);
  assign_el("74",el,2,elec);
  assign_el("75",el,3,elec);

/* write result */


  write_elec(filename);
 }


