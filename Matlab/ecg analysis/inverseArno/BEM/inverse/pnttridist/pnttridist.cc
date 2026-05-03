
  /*************************************************************************/
  /*									   */
  /*	pnttridist: program to compute distance between a point 	   */
  /*	and a triangluated surface.					   */
  /*									   */
  /*									   */
  /* 	20-10-00				Thom Oostendorp		   */
  /*						Medical Physics		   */
  /*						University of Nijmegen	   */
  /*									   */
  /*************************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "trilib.h"
/*
#include "thom.h"
*/

#define DELTA 0.01

void outofmem(void)
 {
  printf("Out of meomry\n");
  exit(1);
 }

/* pntontotri projects r onto triangulated surface */

int pntontotri(double *r, int npnt, int ndhk, double (*pnt)[3], int (*dhk)[3])
 {
  int	idhk, k, mindhk=0;
  double	dist, mindist=0, rp[3], rpmin[3], lab, mu;
  
  for (idhk=0; idhk<ndhk; idhk++)
   {
    dist=ptridist(r, pnt[dhk[idhk][0]], pnt[dhk[idhk][1]], 
		pnt[dhk[idhk][2]], rp, &lab, &mu);
    if (idhk==0 || dist<mindist)
     {
      mindhk=idhk;
      mindist=dist;
      for (k=0; k<3; k++)
        rpmin[k]=rp[k];
     }
   }
  for (k=0; k<3; k++)
    r[k]=rpmin[k];
  return mindhk;
 }


void print_help(void)
 {
  puts("");
  puts("  usage: pnttridist [-options] <tri_file> <pnt_inp_file>");
  puts("");
  puts("    options: -v   be less verbose");
  puts("             -w   be even less verbose");
  puts("             -h   print this help text");
  puts("");
  puts("  The distances between the points in <pnt_inp_file>");
  puts("  and the surface described by <tri_file> are computed ");
  puts("");
  puts("  The output is written to standard output.");
  puts("  If verbose, all individual distances are reported;");
  puts("  in any case the statistics of the distances are reported.");
  puts("");
  exit(1);
 }

bool lessVerbose=false;

void get_arg(int argc, char *argv[], char dhkname[], char pntinname[])
 {
  int	i=1;
  double f;
  
  while (i<argc)
   {
    if (argv[i][0]=='-')
     {
      switch(tolower(argv[i][1]))
       {
case 'v':
	set_verbose(false);
	break;
case 'w':
	set_verbose(false);
	lessVerbose=true;
	break;
default:
	print_help();
       }
     }
    else if (dhkname[0]=='\0')
      strcpy(dhkname, argv[i]);
    else if (pntinname[0]=='\0')
      strcpy(pntinname, argv[i]);
    else
      print_help();
    i++;
   }
   if (pntinname[0]=='\0')
     print_help();
 }


int main(int argc, char *argv[])
 {
  int 	i, k, l, idhk, npnt, npnt1, ndhk, (*dhk)[3], ontr, inside, first=1;
  int	nin=0, nout=0, ntouch=0;
  double (*pnt)[3], (*pnt1)[3], (*pnt2)[3], d;
  double total, r[3][3], mindist=0, maxdist=0, meandist=0, rmsdist=0;
  char	dhkname[80]="\0";
  char	pntinname[80]="\0";
  char	s[80];


  get_arg(argc, argv, dhkname, pntinname);

  if (!dhkinp(dhkname, &npnt, &ndhk, &pnt, &dhk))
    exit(1);
  if (!pntinp(pntinname, &npnt1, &pnt1))
    exit(1);

  pnt2=(double(*)[3])calloc(npnt1, 3*sizeof(double));
  if (pnt2==NULL)
    outofmem();
  for (i=0; i<npnt1; i++)
    for (k=0; k<3; k++)
      pnt2[i][k]=pnt1[i][k];

  if (is_verbose()) printf("\n");
  if (lessVerbose)
    printf ("%s %s ", pntinname, dhkname);
  else
    printf ("Distances between points of %s and surface %s:\n",
	    pntinname, dhkname);
  if (is_verbose()) printf("\n");
  for (i=0; i<npnt1; i++)
   {
    pntontotri(pnt2[i], npnt, ndhk, pnt, dhk);
    d=vecdist(pnt1[i], pnt2[i]);
    meandist += d;
    rmsdist += d*d;
    if (i==0 || d<mindist)
      mindist=d;
    if (i==0 || d>maxdist)
      maxdist=d;
    if (is_verbose())
     {
      total=0;
      for (idhk=0; idhk<ndhk; idhk++)
       {
	for (k=0; k<3; k++)
	  for (l=0; l<3; l++)
	    r[k][l]=pnt[dhk[idhk][k]][l]-pnt1[i][l];
	total += rhoek(r[0], r[1], r[2], &ontr);
	if (ontr)
	 {
	  inside=2;
	  goto finished;
	 }
       }
      if (total>=-DELTA && total <= DELTA)
	inside=0;
      else if (total>=4*M_PI-DELTA && total <= 4*M_PI+DELTA)
	inside=1;
      else
       {
	inside=3;
	if (first && is_verbose())
	 {
	  printf("\nSurface %s is not closed!\n\n", dhkname);
	  first=0;
	 }
       }
finished:
      switch (inside)
       {
case 0:
        strcpy(s, "(outside)");
	nout++;
	break;
case 1:
        strcpy(s, "(inside)");
	nin++;
	break;
case 2:
        strcpy(s, "(touches)");
	ntouch++;
	break;
case 3:
        strcpy(s, "\0");
	break;
       }
      printf("%4d: %8.4f %s\n", i+1, d, s);
     }
   }
  if (lessVerbose)
   {
    printf("%8.4f %8.4f %8.4f %8.4f\n",
	   mindist, maxdist, meandist/npnt1, sqrt(rmsdist/npnt1));
   }
  else
   {
    if (is_verbose())
      printf("\n%4d points inside; %4d points outside; %4d points touch\n",
	     nin, nout, ntouch);
    printf("min dist: %8.4f; max dist: %8.4f; mean dist: %8.4f; RMS dist: %8.4f\n",
	   mindist, maxdist, meandist/npnt1, sqrt(rmsdist/npnt1));
    if (is_verbose()) printf("\n");
   }
  return 0;
 }

