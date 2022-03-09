#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include "planets.h"


/**-----------------Para cargar desde el inicio sistema en coordenadas esféricas (opc 0)**/
#define PI 4.*atan(1.)
#define PI2 pow(PI,2)
#define G 1.0
#define M 1.0
/**----------------Para sistema solar, tiempo en años, distancia en UA----------------**/
#define UAM 1.0/1.495978707e11
#define YRSE 365.25*24*3600.0

//#define GM 1.32712440018e20*pow(UAM,3)*pow(YRSE,2) Es lo mismo
#define GM 39.4769264142519
/**---------------Opción  0 anterior otras unidades (No se usan)---------------------------------------------------------------------**/			
#define DT 0.0001
#define MP 1e-9
#define MT 0.001
#define RP 0.00153172431509222 //Asumiendo misma densidad que titan y masa MP

/**--------Masas de los planetas exteriores relativas al sol (mass in simulation units)--------------**/
#define MJ 9.54791938e-4 
#define MS 2.85885980e-4
#define MU 4.36624404e-5
#define MN 5.15138902e-5
/**--------------------------------------------------------------------------------------------------**/



/** Allocate and free square matrices **/

double **square_matrix(int n) {
  double **mat;
  int i;
  mat = malloc(n * sizeof(double *));
  for(i = 0; i < n; i++)
    mat[i] = malloc(n * sizeof(double));
  return mat;
}
void free_square_matrix(double **mat, int n) {
  int i;
  for(i = 0; i < n ; i++)
    free(mat[i]);
  free(mat);
}

double **nonsquare_matrix(int n,int m) {
  double **mat;
  int i;
  mat = malloc(m * sizeof(double *));
  for(i = 0; i < m; i++)
    mat[i] = malloc(n * sizeof(double));
  return mat;
}



/** Funcion que realiza el producto de el vector por la matriz**/

double *mat_vect(double *v, double **mat,int n) {
  double *vt;
  int i, j;
  vt = malloc(n*sizeof(double));
  if(i == j)
    for(i = 0; i < n; i++)
      vt[i] = 0;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      vt[i]+=mat[i][j]*v[j];
  return(vt);
}
/** Funcion que realiza el producto de un matriz por otra **/ double

**mat_mat(double **mat1, double **mat2,int n) {
  double **mat3;
  int i, j, m;
  mat3=square_matrix(n);
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      mat3[i][j] = 0;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      for(m = 0; m < n; m++)
	mat3[i][j]+=mat1[i][m]*mat2[m][j];
  return(mat3);
}

/**esta funcion me muestra la matriz**/
void print_matrix(double **mat, int n) {
  int i,j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      printf("%1.20f\n", mat[i][j]);
    printf("\n");
  }
}
void print_vect(double *vec, int n, char *fname) {
  int i,j,fdesc;
  char filename[FNAMESIZE + 6] = "datos/";
  FILE *arch1;
  strcat(filename, fname);
  if((arch1=fopen(filename, "w")) == NULL){
    printf("No se pudo abrir/crear el archivo %s\n",filename);
  }
  
  else
    printf(" - Datos guardados en %s\n", filename);
  for(i = 0; i < n; i++)
    fprintf(arch1,"%d\t%f\n", i,vec[i]);
  printf("\n");
  close(filename);
}




/**----------------------------------------------------------------------------------------------------------------------------**/
/**              Carga los elementos orbitales de los planetas y partículas desde el archivo "sist_solar.dat" **/
/**------------------------------------------------------------------------------------------------------------------------------**/
void elem_orbitales(SIS *sis, PAR *par){
  FILE *arch;
  int i,exp,j,n,nplan;
  double t;
  
  double aa[7];
  /* Abrir elementos  orbitales */
  if((arch=fopen("sist_solar.dat", "r")) == NULL){
    printf("No se puede leer el archivo/%s\n","sist_solar.dat");
  }
  else
    printf(" - Datos de elementos orbitales cargados\n");
  while(fscanf(arch,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",&i,&aa[0],&aa[1],&aa[2],&aa[3],&aa[4],&aa[5],&aa[6],&exp)!=EOF){
  	
    aa[6]=aa[6]/(pow(10.,exp));
	printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%1.9f\n",i,aa[0],aa[1],aa[2],aa[3],aa[4],aa[5],aa[6]);
    sis->x[i] = aa[0];
    sis->y[i] = aa[1];
    sis->z[i] = aa[2]*2.*PI/360.;	
    sis->vx[i] = aa[3]*2.*PI/360.;
    sis->vy[i] = aa[4]*2.*PI/360.;
    sis->vz[i] = aa[5]*2.*PI/360.;
    sis->m[i] = aa[6];
	
 }
  fclose(arch);

}

/** -----------------------------------------------------------------------**/ 
int guardar_datos(SIS *r, PAR *par, char *fname) {
  char filename[FNAMESIZE + 6] = "datos/";
  int fdesc;
  int n = (par->n)+1;
  double t = par->time;
  
 
  strcat(filename, fname);
  
  fdesc = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (fdesc < 0) {
    printf("*** Cannot open %s.\n", filename);
    return 0;
  }	
  

  write(fdesc, par, sizeof(PAR));
  write(fdesc, r->x, n * sizeof(double));	
  write(fdesc, r->y, n * sizeof(double));	
  write(fdesc, r->z, n * sizeof(double));
  write(fdesc, r->vx, n * sizeof(double));	
  write(fdesc, r->vy, n * sizeof(double));
  write(fdesc, r->vz, n * sizeof(double));
  write(fdesc, r->m, n * sizeof(double));
  close(fdesc);
  return 1;
}


void print_sistem(SIS *sis,PAR *par, char *fname,double time) {
  int n,i;
 
  char filename[FNAMESIZE + 6] = "datos/";
  FILE *arch1;	
  n=par->n;
 
  strcat(filename, fname);
  if((arch1=fopen(filename, "w")) == NULL){
    printf("No se pudo abrir/crear el archivo %s\n",filename);
  }
  
  else
    printf(" - Datos guardados en %s\n", filename);
    fprintf(arch1,"%1.3f\n",time);

    for (i=0;i<=n;i++){
        fprintf(arch1,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i, sis->x[i],sis->y[i],sis->z[i],sis->vx[i],sis->vy[i],sis->vz[i],sis->m[i]);
    }
  printf("\n");
  fclose(arch1);
 
}


void inicializar_sistema(SIS *sis, PAR *par, int opc, char *fname){
  
  double ecc,x0,dx,**rot,**rxq,*vq,*vrect,*rrect,step,nmedio;
  double r;
  char filename[FNAMESIZE];
  double theta,phi,dot_theta,dot_phi, delta,accu,*q,q1,q2,*dqdt;  
  char *p;
  int fdesc, count;
  int n,hn,t,i,j,k;
  //double **fitr,**fitv;
 
  n=(par->n)+1;
  q=malloc(3 * sizeof(double));
  dqdt=malloc(3 * sizeof(double));
  vq= malloc(6 * sizeof(double));
  vrect = malloc(3 * sizeof(double));
  rrect= malloc(3 * sizeof(double));
  x0=0.6;
  dx=1.8/n;
 step=1.e-6;
 accu=1.e-14;
 hn=(int) sqrt(n);
 rot=square_matrix(3);
 rxq=square_matrix(3);

init_ran(10);
 if (opc==0){
	i=1;
	sis->x[0]=0.0;
   	sis->y[0]=0.0;
   	sis->z[0]=0.0;
   	sis->vx[0]=0.0;
   	sis->vy[0]=0.0;
   	sis->vz[0]=0.0; 
   	sis->m[0]=1.0;	
  for (j=1; j<=hn;j++){
	for (k=1; k<=hn;k++){
	//printf("%d\t",i);
	theta=(PI/hn)*(j-1);
	phi=((2.*PI)/hn)*(k-1);
	r=1.0+(dran()*1.e-9);
	if (r<0.0){
		r=-r;
		
	}
	
	//printf("%f\n",r);
	sis->x[i]=(double) r*sin(theta)*cos(phi);
        sis->y[i]=(double) r*sin(theta)*sin(phi);
	sis->z[i]=(double) r*cos(theta);
	vrect[0]=0.0;
	dot_theta=sqrt(GM/(2.*(1.+sin(theta))*pow(r,3)));
	vrect[1]=(double) r*dot_theta;
	vrect[2]=(double) r*sin(theta)*dot_theta;
	rot[0][0]=sin(theta)*cos(phi);
	rot[0][1]=cos(theta)*cos(phi);
	rot[0][2]=-sin(phi);
	rot[1][0]=sin(theta)*sin(phi);
        rot[1][1]=cos(theta)*sin(phi);
	rot[1][2]=cos(phi);
	rot[2][0]=cos(theta);
	rot[2][1]=-sin(theta);
	rot[2][2]=0.0;
	vrect=mat_vect(vq, rot,3);
	 sis->vx[i]=vrect[0];
         sis->vy[i]=vrect[1];
	 sis->vz[i]=vrect[2];
 	 sis->m[i]=MP;
	i++;
 	}
	}
	}
	else if(opc==1){

	
	/**Uso variable sis (x,y,z,vx,vy,vz) por comodidad pero estoy leyendo (a,e,i,w,W,M)**/

	elem_orbitales(sis,par);//Leo elementos orbitales de sist_solar.dat 
	sis->x[0]=0.0;
   	sis->y[0]=0.0;
   	sis->z[0]=0.0;
   	sis->vx[0]=0.0;//0.00202619;
   	sis->vy[0]=0.0;//-0.00462527;
   	sis->vz[0]=0.0;//1.2715e-05; 
	sis->m[0]=1.0;
		
	for (i=1;i<=par->n;i++){
		vq[0]=sis->x[i];//a
		vq[1]=sis->y[i];//e
		vq[2]=sis->z[i];//i
		vq[3]=sis->vx[i];//w   
		vq[4]=sis->vy[i];//W
		vq[5]=sis->vz[i];//M
		printf("%d\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.9f\n",i, sis->x[i],sis->y[i],sis->z[i],sis->vx[i],sis->vy[i],sis->vz[i],sis->m[i]);

	/**------------------Matriz de cambio de base Rxq para pasar de coordenadas q=(q1,q2,0) -----> r=(x,y,z)-----**/
	/**-------------------q1=a*(cos(E)-e); q2=a*sqrt(1-e²)*sin(E) -----------------------------------------------**/

	rxq[0][0]=(cos(vq[4])*cos(vq[3]))-(sin(vq[4])*cos(vq[2])*sin(vq[3]));
	rxq[0][1]=(-cos(vq[4])*sin(vq[3]))-(sin(vq[4])*cos(vq[2])*cos(vq[3]));
	rxq[0][2]=sin(vq[4])*sin(vq[2]);
	rxq[1][0]=(sin(vq[4])*cos(vq[3]))+(cos(vq[4])*cos(vq[2])*sin(vq[3]));
        rxq[1][1]=(-sin(vq[4])*sin(vq[3]))+(cos(vq[4])*cos(vq[3])*cos(vq[2]));
	rxq[1][2]=-cos(vq[4])*sin(vq[2]);
	rxq[2][0]=sin(vq[2])*sin(vq[3]);
	rxq[2][1]=sin(vq[2])*cos(vq[3]); 
	rxq[2][2]=cos(vq[2]);

	/**------------------Para esto hay que resolver antes la ecuación de Kepler: --------------------------------**/
	/**--------------------------------   E-(e*sin(E))= n(t-t0) = M (vq[5]!!!)  ---------------------------------**/
	ecc=0.0;
	delta=(ecc-(vq[1]*sin(ecc)))-vq[5];
	while(delta < 0.0){
	ecc+=step;
	delta=(ecc-(vq[1]*sin(ecc)))-vq[5];	
	}
	printf("%1.10f\n",ecc/(2.*PI/360.));
	q1=vq[0]*(cos(ecc)-vq[1]);
	q2=vq[0]*sqrt(1.-pow(vq[1],2))*sin(ecc);
	q[0]=q1;
	q[1]=q2;
	q[2]=0.0;
	
	/**-------------------------Ahora sí -------------------------------------------------**/
	/**-------------------------Coordenadas rectangulares posición-----------------------------------**/

	rrect=mat_vect(q, rxq,3);
	sis->x[i]=rrect[0];
	sis->y[i]=rrect[1];
	sis->z[i]=rrect[2];


	/**-----------------------Coordenadas de velocidad----------------------------------**/

	nmedio=sqrt(GM*(sis->m[0]+sis->m[i])/pow(vq[0],3)); // M =n (t-t0), n = sqrt(G*(m0+m1))a**(-3/2) 3ra Kepler
	printf("%d\t%f\n",i,nmedio);
	/** --------------------------dq1/dt = -n*a*sin(E)/(1-e*cos(E))-----------------------**/
	/** --------------------------dq2/dt = sqrt(1-e**2)*cos(E)/(1-e*cos(E)) --------------**/	

	dqdt[0]= -nmedio*vq[0]*sin(ecc)/(1.-(vq[1]*cos(ecc)));
	dqdt[1]=nmedio*vq[0]*sqrt(1.-pow(vq[1],2))*cos(ecc)/(1.-(vq[1]*cos(ecc)));
	dqdt[2]=0.0;

	vrect=mat_vect(dqdt, rxq,3);
		
	sis->vx[i]=vrect[0];
	sis->vy[i]=vrect[1];
	sis->vz[i]=vrect[2];
		
	printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,sis->x[i],sis->y[i],sis->z[i],sis->vx[i],sis->vy[i],sis->vz[i]);
	}	
	}	
	
else if (opc==-1){
     
    fdesc = open(fname, O_RDONLY, 0644);
    if (fdesc < 0) {
      printf("--- Read %s not successful.\n", fname);
      return 0;
    }
	
    else{

      read(fdesc, par, sizeof(PAR)); 
      
       n = (par->n)+1;
	
	sis->x=malloc((n)*sizeof(double));
	sis->y=malloc((n)*sizeof(double));
 	sis->z=malloc((n)*sizeof(double));
	sis->vx=malloc((n)*sizeof(double));
  	sis->vy=malloc((n)*sizeof(double));
	sis->vz=malloc((n)*sizeof(double));
  	sis->m=malloc((n)*sizeof(double));
    
      	read(fdesc, sis->x,(n)*sizeof(double));	
	read(fdesc, sis->y,(n)*sizeof(double));    
	read(fdesc, sis->z,(n)*sizeof(double));    
        read(fdesc, sis->vx,(n)*sizeof(double));
	read(fdesc, sis->vy,(n)*sizeof(double));
	read(fdesc, sis->vz,(n)*sizeof(double));
        read(fdesc, sis->m,(n)*sizeof(double));
	close(fdesc);
      }
 }

}


void verlet_grav_3d(SIS *sis,PAR *par){
  int t, n, nplan,k,l;
  double tstep,rp,rm;
  double time;  
  
 // double *xx,*yy,*zz,*vxx,*vyy,*vzz,*mm;
  //double *fx,*fy,*fz,*fvx,*fvy,*fvz,*mp;
  n=(par->n)+1;
	
  nplan=par->np;
  time=par->time;
  double x0,dx,r;
  double fGM, rij;
  int i,j;
  double *x,*y,*z,*ax,*ay,*az,*ax_,*ay_,*az_;
  x=malloc(n*sizeof(double));
  y=malloc(n*sizeof(double));
  z=malloc(n*sizeof(double));
  ax=malloc(n*sizeof(double));
  ay=malloc(n*sizeof(double));
  az=malloc(n*sizeof(double));
  ax_=malloc(n*sizeof(double));
  ay_=malloc(n*sizeof(double));
  az_=malloc(n*sizeof(double));
 
  x0=0.6;
  dx=1.8/(n-3);
  tstep=DT;


/**------------------------------Acá empezaría si la ultima opc (tnp) es 0 --------------------------------------------------- **/
  
	

	for (i=0; i<n;i++){

  
        ax[i]=0.0;
        ay[i]=0.0;
        az[i]=0.0;

	
	for (j=0;(j<=nplan);j++){
	if(i==j){
	ax[i]+=0.0;
	ay[i]+=0.0;
  	az[i]+=0.0;
	}
	else{
	rij=sqrt(pow(sis->x[i]-sis->x[j],2)+pow(sis->y[i]-sis->y[j],2)+pow(sis->z[i]-sis->z[j],2));

  	fGM=GM*(sis->m[j])/pow(rij,2);
        ax[i]+=fGM*(sis->x[j]-sis->x[i])/rij;
	ay[i]+=fGM*(sis->y[j]-sis->y[i])/rij;
  	az[i]+=fGM*(sis->z[j]-sis->z[i])/rij;
	}
	}
        x[i]=sis->x[i]+(sis->vx[i]*tstep)+(ax[i]*(pow(tstep,2))/2.);
        y[i]=sis->y[i]+(sis->vy[i]*tstep)+(ay[i]*(pow(tstep,2))/2.);
        z[i]=sis->z[i]+(sis->vz[i]*tstep)+(az[i]*(pow(tstep,2))/2.);
	}


	for (i=0; i<n;i++){
        ax_[i]=0.0;
        ay_[i]=0.0;
        az_[i]=0.0; 
    	
    	for (j=0;(j<=nplan);j++){
	if(i==j){
	ax[i]+=0.0;
	ay[i]+=0.0;
  	az[i]+=0.0;
	}
	else{
        rij=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
        
	fGM=GM*(sis->m[j])/pow(rij,2);    
            
	  ax_[i]+=fGM*(x[j]-x[i])/rij;
	  ay_[i]+=fGM*(y[j]-y[i])/rij;
	  az_[i]+=fGM*(z[j]-z[i])/rij;
	}
	
	}
	
    }
    for (i=0; i<n;i++){
	sis->x[i]=x[i];
	sis->y[i]=y[i];
	sis->z[i]=z[i];
	
  	sis->vx[i]=sis->vx[i]+((ax[i]+ax_[i])*(tstep/2.)); 
  	sis->vy[i]=sis->vy[i]+((ay[i]+ay_[i])*(tstep/2.));
        sis->vz[i]=sis->vz[i]+((az[i]+az_[i])*(tstep/2.));
    }
	
free(x);
free(y);
free(z);
free(ax);
free(ay);
free(az);
free(ax_);
free(ay_);
free(az_);
    
}

/*************************************************************************************************************************************************/

/** Programa principal --- Leo parámetros desde la ventana de comandos y guardo configuración inicial**/

/**************************************************************************************************************************************************/

int main(int argc, char *argv[]) {
  SIS sis;

  PAR par;
  int nsteps,i,l,j,k,jmax;
  int opc,nplan;
  double count1,count2,time,frec1,frec2;
  char *filename;//[FNAMESIZE];  
  char fname[FNAMESIZE];
 

  if (argc < 7) {
    printf("Uso: %s <N particles> <Tmax> <option(0,1,-1)> <file> <time_write_config> <nr planets> \n", argv[0]);
    exit(1);
}
  // Lee los parametros requeridos desde la linea de comandos.
  par.n = atoi(argv[1]);
  par.tmax = atof(argv[2]);
  opc=atoi(argv[3]);
  filename=argv[4]; 
  frec1=atof(argv[5]);
  par.np=atoi(argv[6]);


  
  
  	
  
  sis.x=malloc(((par.n)+1)*sizeof(double));
  sis.y=malloc(((par.n)+1)*sizeof(double));
  sis.z=malloc(((par.n)+1)*sizeof(double));

  sis.vx=malloc(((par.n)+1)*sizeof(double));
  sis.vy=malloc(((par.n)+1)*sizeof(double));
  sis.vz=malloc(((par.n)+1)*sizeof(double));
  sis.m=malloc(((par.n)+1)*sizeof(double));
  
  printf("\n");
  printf("-----------------------------------------------------------------------------------------------------\n");
  printf(" --- Evolución de sistema gravitatorio --> Masa central y N particulas \n");
  printf("-----------------------------------------------------------------------------------------------------\n");
  printf("\n");
  printf("-----------------------------\n");
  printf(" Parametros ingresados \n");
  printf("-----------------------------\n");
  printf(" Nparticulas = %s\n", argv[1]);
  printf(" Tiempo de integración = %s\n", argv[2]);
  printf("------------------------------\n");
  printf("\n");

  /**----------------------------Inicializo sistema-----------------------------------------**/
  inicializar_sistema(&sis,&par,opc,filename);
  /**-------------------------------------------------------------------------------------------------------**/
  par.tmax = atof(argv[2]);
 // par.np=atoi(argv[7]);
 // par.tnp=atoi(argv[8]);

  printf("\n");
 
  time=par.time;
 
count1=0.0;


if(opc!=-1){
par.time=0.0;
time=par.time;
printf("%d\t%f\t%f\%s\n",par.n,par.tmax,frec1,filename);
}
else{
time=par.time;
}

while(time <= par.tmax){
		
		verlet_grav_3d(&sis,&par);
	        
               	if (count1 >= frec1){
		
	 	sprintf(fname, "%d_%lf_%lf", par.n,par.tmax,time);
		printf("%s\n",fname);
		guardar_datos(&sis,&par,fname);
		count1=0.0;
		}
		count1+=DT;
		time+=DT;
		par.time=time;
	  
}
	
free(sis.x);
free(sis.y);
free(sis.z);
free(sis.vx);
free(sis.vy);
free(sis.vz);

  exit(EXIT_SUCCESS);


}
