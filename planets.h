// Maximum file name size:
#define FNAMESIZE 64



typedef struct PAR {
	int n;
	int np;
	double tmax;
	double time;
	int tnp;
   
} PAR;

/**typedef struct FIT {
  
  double **jr, **jv;
  double **sr, **sv;
  double **ur,**uv;
  double **nr,**nv;

  
} FIT;
**/

typedef struct SIS {
  double *x, *y, *z;
  double *vx, *vy, *vz;
  double *m;

} SIS;


