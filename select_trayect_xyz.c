#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include "planets.h"

void result(PAR *par,SIS *sis)
{
 int n,i; 
 double r,ec;
 FILE *arch1,*arch2;
 char fname1[FNAMESIZE],fname2[FNAMESIZE]; 
 n=par->n;
 i=0;
 //sprintf(fname1, "radial.dat");
 //sprintf(fname2, "sistem.dat");
/**
if((arch1=fopen(fname1, "w")) == NULL){
   printf("No se pudo abrir/crear el archivo %s\n",fname1);
  }
  
if((arch2=fopen(fname2, "w")) == NULL){
    printf("No se pudo abrir/crear el archivo %s\n",fname2);
  }**/
    //for (i=0;i<=n;i++){
	r=sqrt(pow(sis->x[i]-sis->x[0],2)+pow(sis->y[i]-sis->y[0],2)+pow(sis->z[i]-sis->z[0],2));
	ec=(pow(sis->vx[i],2)+pow(sis->vy[i],2)+pow(sis->vz[i],2))*sis->m[i]/2.;
	//printf("%lf\n",);
	//printf("%lf\t%lf\t%lf\n ", par->time,r,ec);	 
	printf(" %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n ",par->time,sis->x[i]-sis->x[0],sis->y[i]-sis->y[0],sis->z[i]-sis->z[0],sis->vx[i]-sis->vx[0],sis->vy[i]-sis->vy[0],sis->vz[i]-sis->vz[0]),r,ec;
	  //printf(" %d\t%1.30f\t%1.30f\n ",i,x[i],y[i]);

//}
    //printf(" - Datos guardados en %s\n", fname1);
    //printf(" - Datos guardados en %s\n", fname2);
//	printf("%f\t%1.3f\t%1.3f\n",time, x[part],y[part]);//,sis->vx[part],sis->vy[part]);

  //printf("\n");
  //fclose(arch1);
  //fclose(arch2);	
 	
    return;
}
main()
{
  char filename[FNAMESIZE + 5];
  int fdesc;
  int i = 0;
  int j = 0;
  PAR par;
  SIS sis;
 
  char *p;
  int n;
  double time;
  //int part; 

  //part=atoi(argv[2]);
  while(fgets(filename, FNAMESIZE + 5, stdin)){
    j = 0;
    p = strchr(filename, '\n');
   *p = '\0';
 	
    fdesc = open(filename, O_RDONLY, 0644);

    if (fdesc < 0) {
      printf("--- Read %s not successful.\n", filename);
      return 0;
    }
    else{
      
      //read(fdesc,time,sizeof(double));
      read(fdesc, &par, sizeof(PAR)); 
      n = (par.n)+1; 
  	sis.x=malloc((n)*sizeof(double));
  	sis.y=malloc((n)*sizeof(double));
	sis.z=malloc((n)*sizeof(double));
 	 sis.vx=malloc((n)*sizeof(double));
	  sis.vy=malloc((n)*sizeof(double));
 	sis.vz=malloc((n)*sizeof(double));
  	sis.m=malloc((n)*sizeof(double));
 
      read(fdesc, sis.x, n*sizeof(double)); 
      read(fdesc, sis.y, n*sizeof(double)); 
	read(fdesc, sis.z, n*sizeof(double)); 
      read(fdesc, sis.vx, n*sizeof(double)); 
      read(fdesc, sis.vy, n*sizeof(double));     
	read(fdesc, sis.vz, n*sizeof(double));     
      read(fdesc, sis.m, n*sizeof(double)); 


//printf("%s",p);
       result(&par,&sis);
	
      
      close(fdesc);
    }
  }
 
  
 
  return 0;
    
}
  
