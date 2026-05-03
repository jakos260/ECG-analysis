/* file matg.h */

struct L* floggt1 (void);

#define maxpd 30
#define nnbmax 50  /*** be careful !!! ***/
#define maxt 40
#define nsmax 5
#define npmax 800
#define ntmax 1575
#define nphmax 300
#define nthmax 600
#define nptmax 300
#define nttmax 600
#define nelmax 64
#define natmax 151
#define dipstr 40
#define VDIP 40.0
#define rhtol 0.001
#define xtol 0.0001
#define small 0.053674316E-06
#define big 0.104857600E06
#define pi 3.141592653589793


# define max(a,b) 		(a<b ? b : a)
# define min(a,b) 		(a>b ? b : a)
static float maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define MIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) < (maxarg2) ?\
	(maxarg1) : (maxarg2))


typedef struct S
	     {
	       char   qgeonm [10];
	       int    npt ;	     /* number of points of surface [i] */
	       int    ifstpt ;    /* first point of surface [i] */
	       int    ilstpt ;    /* last point of surface [i] */
	       int    ntr;
	       int    ifsttr;
	       int    ilsttr;
	     } surface;

typedef struct H
	      {
	       char   qhrtnm [20];   /* name of the heart */
	       int    npntsh;        /* number of points of heart */
	       int    ntrsh;         /* number of triangles ... */
	      } heart;

typedef struct L
	      {
	       char   qcase [10];    /* nazov pripadu */
	       char   qdate [15];
	       char   qtime [15];
	       int    ns;            /* pocet povrchov */
	       int    npnts;         /* total number of points */
	       int    ntrs;          /* total number of triangles */
	       struct S  povrch [10];
	       struct H  srdce;
	       int    icont [10][10];  /* vzajomne ulozenie povrchov */
	       float  sigpls [10];
	       float  sigmns [10];
	       float  gamma [10][10];
	      } zadanie;

 typedef struct
	  {
	   float x;
	   float y;
	   float z;
	   } suradnice;

 typedef struct Q
	  {
	   int    pd; /*pocet dipolov*/
	   int  dx [maxpd];
	   int  dy [maxpd];
	   int  dz [maxpd];
	   int  stredx;
	   int  stredy;
	   int  stredz;
	   int  polomer;
	  } srd;

