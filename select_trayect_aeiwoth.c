#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include "planets.h"

#define GM 39.4769264142519
#define PI 4.*atan(1.)
/** Funcion producto interno **/
double vect_vect(double *v, double *u,int n) {
  double vu;
  int i, j;
  //vu = malloc(n*sizeof(double));
  vu=0.0;	
    for(i = 0; i < n; i++)
    	vu+=v[i]*u[i];
  return(vu);
}
double *vectxvect(double *v, double *u,int n) {
  double *vu;
  int i, j;
  vu = malloc(n*sizeof(double));
	

    	vu[0]=v[1]*u[2]-v[2]*u[1];
    	vu[1]=-(v[0]*u[2]-v[2]*u[0]);    	
	vu[2]=v[0]*u[1]-v[1]*u[0];

  return(vu);
}
void result(PAR *par,SIS *sis)
{
 int n,i,j; 
 double theta,r,ec,a,e,inc,w,o,l,g,omm,v,mn,mh,mu,ham,sigma,*h,*vecr,*vece,ne,*veck,*vecn,*vecv,*angl,*pangl,*c,*vecc,*vecp;
 FILE *arch1,*arch2;
 char fname1[FNAMESIZE],fname2[FNAMESIZE]; 
 n=(par->n);
 i=1;
 vecr = malloc(3*sizeof(double));
 vecv = malloc(3*sizeof(double));
 h=malloc(3*sizeof(double));
 vece = malloc(3*sizeof(double));
 vecc = malloc(3*sizeof(double));
 vecn = malloc(3*sizeof(double));
 veck = malloc(3*sizeof(double));
 angl=malloc(3*sizeof(double));
 pangl=malloc(3*sizeof(double));
 c = malloc(3*sizeof(double));
 vecp = malloc(3*sizeof(double));
 //sprintf(fname1, "radial.dat");
 //sprintf(fname2, "sistem.dat");
 
    //for (i=1;i<=n;i++){
/**	
	 sprintf(fname1, "%d_rect.dat",i);
	 sprintf(fname2, "%d_orb.dat",i);
	if((arch1=fopen(fname1, "w")) == NULL){
   	printf("No se pudo abrir/crear el archivo %s\n",fname1);
 	 	}
  
	 if((arch2=fopen(fname2, "w")) == NULL){
	    printf("No se pudo abrir/crear el archivo %s\n",fname2);
	  }
   **/
		vecr[0]=sis->x[i]-sis->x[0];
		vecr[1]=sis->y[i]-sis->y[0];
		vecr[2]=sis->z[i]-sis->z[0];
		r=sqrt(vect_vect(vecr,vecr,3));
		vecv[0]=sis->vx[i]-sis->vx[0];
		vecv[1]=sis->vy[i]-sis->vy[0];
		vecv[2]=sis->vz[i]-sis->vz[0];
		veck[0]=0.0;
		veck[1]=0.0;
		veck[2]=1.0;

		
		//for (j=0;j<3;j++){
		//	vecp[j]=(sis->m[i])*vecv[j];
		//	vecc[j]=(mu*(sis->m[i])*vecr[j])/r;
		//	}
		v=sqrt(vect_vect(vecv,vecv,3));
//	r=sqrt(pow(sis->x[i]-sis->x[0],2)+pow(sis->y[i]-sis->y[0],2)+pow(sis->z[i]-sis->z[0],2));
//	v=sqrt(pow(sis->vx[i],2)+pow(sis->vy[i],2)+pow(sis->vz[i],2));
		ec=(pow(v,2))/2.;
		ham=ec-(GM*(sis->m[0]+sis->m[i])/r);
		mu=GM*(sis->m[0]+sis->m[i]);
		
		a=-mu/(2.*ham);
		h=vectxvect(vecr,vecv,3);
		mh=sqrt(vect_vect(h,h,3));
		angl=vectxvect(vecr,vecv,3);
		
		e=sqrt(1.-(pow(mh,2)/(mu*a)));
		vecc=vectxvect(vecv,h,3);
		for(i=0;i<3;i++){
			vecc[i]=vecc[i]-(mu*vecr[i]/r);
			vece[i]=vecc[i]/mu;		
		}
		inc=acos(h[2]/mh)*360./(2.*PI);
		vecn=vectxvect(veck,h,3);
		mn=sqrt(vect_vect(vecn,vecn,3));
		if(vecn[1]>=0.0)		
		omm=acos(vecn[0]/mn)*360./(2.*PI);
		else
		omm=360.-(acos(vecn[0]/mn)*360./(2.*PI));
		ne=vect_vect(vecn,vece,3);
		if(vece[2]<0.0){
			w=360.-(acos(ne/(mn*e))*360./(2.*PI));
		}		
		else{
			w=(acos(ne/(mn*e))*360./(2.*PI));
		}
		
		if(vect_vect(vecr,vecv,3)>=0.0){
			theta=acos(vect_vect(vece,vecr,3)/(e*r))*360./(2.*PI);
		}
		else{
			theta=360.-(acos(vect_vect(vece,vecr,3)/(e*r))*360./(2.*PI));
		}
		//fprintf(arch1," %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n ",par->time,sis->x[i]-sis->x[0],sis->y[i]-sis->y[0],sis->z[i]-sis->z[0],sis->vx[i]-sis->vx[0],sis->vy[i]-sis->vy[0],sis->vz[i]-sis->vz[0]);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",par->time,a,e,inc,w,omm,theta);  
		
	 //fclose(arch1);
  	 //fclose(arch2);
	//}
//}
    //printf(" - Datos guardados en %s\n", fname1);
    //printf(" - Datos guardados en %s\n", fname2);
//	printf("%f\t%1.3f\t%1.3f\n",time, x[part],y[part]);//,sis->vx[part],sis->vy[part]);

  //printf("\n");
 	
 free(vecr);
 free(vecv);
 free(h);
 free(vece);
 free(vecc);
 free(vecn);
 free(veck);
 free(angl);
 free(pangl);
 free(c); 
 free(vecp);
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
     //return 0;
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
	
      
      
    }
  }
 
  close(fdesc);
 
  return 0;
    
}
  
