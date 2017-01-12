

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
//~ #define NTHREADS 10

struct param{


int N;
int NumPix;
int *resultados; /*No especificar el tamaño e inicializarlo sin más */
double beta;
int contador;
int dix;
int diy;
double deltaPaso;
double deltaX;
double deltaY;
int factordim;
double mlenwid; /*Tamaño de la lente principal*/
int mlendiv; /*Num. de divisiones de la lente principal*/
int nulen; /*Num de lentillas*/
double ff; /*distancia lenslet-detector*/
double fociii; 
double dist;
double anchura;
double angmax;
double *obj;
double deltadis;
double d;
};

//////////////////////////////////////////////
/*Inicializar las funciones*/
double *prop(double delta, double *po1, double *po0);
double *crea_array(int dd);
void guarda_array(double *arra, int dix, int diy, FILE *fuuu);
void trataobj(FILE *fui, int fac, int dix, int diy, int deltax, int deltay); /*El fac2 es solo para este tip de objeto*/
//~ double *crea_matriz(int dx, int dy, int dex, int dey);
 int *crea_array_ent(int dil);
void guarda_resultados(int *arra, int nmpix,  FILE *fiiii);
void guarda_resul(double** array, int NumRayos, FILE *finp, double dist1, double dist2);
FILE *apertura_archivo(void);
void *propagacion(void *threadarg);
double *calcula_mu(int longi, double *coef, FILE *finp);
int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet);
/*double *matrixmul(double aaa, double bbb);*/
/*#define anchura_malla 1
#define deltaY 0.001
#define deltaX 0.001
#define dimx anchura_malla/deltaX
#define dimy anchura_malla/deltaY*/
#define PI 3.141592653
#define INPUT_TXT
#define RAND_ANG
#define SQUARE
#define DOSPI/*La otra opcion es abanico*/
int main(int argc, char **argv)
{
int i, cont, dimx, dimy, factordim, NumAdq,rc, dimtotal, a, b, nth, factorlente;
int array_inp[14];
double array_inp2[16];
double  d, deltaX, deltaY, anchura_malla, deltaang, beta,  angmax, radio, deltax, deltay; /* focii = distancia + anchura_malla/2 ¿? */
//~ pthread_t threads[NTHREADS];
//~ int thread_arg[NTHREADS];
clock_t comienzo, fin;
FILE  *fu1, *finp;
double factorcon=PI/180;
nth=0;
factorlente=10;  /*Numero de detectores por cada minilente*/

#ifdef ABANICO /*Solo se proyecta en un abanico de angulos*/

double angulos[2]={ -45*factorcon, 0*factorcon};
finp=fopen("AngulosBundle.txt", "w");
NumAdq=sizeof(angulos)/sizeof(angulos[0]);
		for(cont=0; cont<NumAdq; cont++)
		{
			fprintf(finp, "%lf\n", angulos[cont]);
		}
fclose(finp);
finp=fopen("NumAngAbanico.txt", "w");
fprintf(finp, "%d\n", NumAdq);
fclose(finp);
#endif



#ifdef INPUT_TXT /*Error aqui */
		fu1=fopen("inputint2.txt", "r");
		fscanf(fu1, "%d", &array_inp[0]); /*NumAdq*/
		fscanf(fu1, "%d", &array_inp[1]); /*N*/
		fscanf(fu1, "%d", &array_inp[2]); /*NumPix*/
		fscanf(fu1, "%d", &array_inp[3]); /*factordim*/
		fscanf(fu1, "%d", &array_inp[4]); /* divisiones de la main lens*/
		fscanf(fu1, "%d", &array_inp[5]); /*num lens*/
		fclose(fu1);

		fu1=fopen("inputdou2.txt", "r");
		fscanf(fu1, "%lf", &array_inp2[0]); /*Anchura*/
		fscanf(fu1, "%lf", &array_inp2[1]); /*deltaX*/
		fscanf(fu1, "%lf", &array_inp2[2]); /*deltaY*/
		fscanf(fu1, "%lf", &array_inp2[3]); /*deltadis*/
		fscanf(fu1, "%lf", &array_inp2[4]); /*d*/
		fscanf(fu1, "%lf", &array_inp2[5]); /*deltapaso*/
		fscanf(fu1, "%lf", &array_inp2[6]); /* Anchura de lente principal*/
		fscanf(fu1, "%lf", &array_inp2[7]); /*f*/
		fscanf(fu1, "%lf", &array_inp2[8]); /*focii*/
		//~ fscanf(fu1, "%lf", &array_inp2[9]); /*distancia*/
		fclose(fu1);

//~ printf("Input enteros:\n");
		//~ for (cont=0; cont< 4; cont++)
		//~ {
			
			//~ printf("array_input[%d]=%d\n",cont, array_inp[cont]);
		//~ }
	
//~ getchar();

//~ printf("Input dobles:\n");
//~ for (cont=0; cont< 6; cont++)
//~ {
	
	//~ printf("array_input2[%d]=%lf\n",cont, array_inp2[cont]);
//~ }

//~ NumAdq= array_inp[0];

#endif



#ifdef DOSPI
NumAdq=array_inp[0]; /*Numero de angulos de adquisición*/
#endif
struct param par[NumAdq];


deltaang=2*PI/NumAdq;
anchura_malla=1; /*En unidades de medida, anchura del objeto*/

srand(time(NULL));
//~ angmax=; 

	for(i=0; i<NumAdq; i++){
			#ifdef INPUT_TXT
			par[i].N=array_inp[1];
			par[i].NumPix=array_inp[2];
			par[i].anchura=array_inp2[0];
			par[i].deltaX=array_inp2[1];
			par[i].deltaY=array_inp2[2];
			par[i].deltadis=array_inp2[3];
			par[i].d=array_inp2[4];
			par[i].factordim=array_inp[3];
			par[i].mlendiv=array_inp[4];
			par[i].nulen=array_inp[5];
			par[i].deltaPaso=array_inp2[5];
			par[i].mlenwid=array_inp2[6];
			par[i].ff=array_inp2[7];
			par[i].fociii=array_inp2[8];
			//~ par[i].dist=array_inp2[9];
			#endif
		
			#ifdef ARGS /*  Faltarán parte de los argumentos */
			par[i].N=atoi(argv[2]);
			par[i].NumPix=atoi(argv[3]);
			par[i].anchura=(argv[4]);
			par[i].deltaX=(argv[5]);
			par[i].deltaY=(argv[6]);
			par[i].deltadis=(argv[7]);
			par[i].d=(argv[8]);
			par[i].factordim=(argv[9]);
			par[i].deltaPaso=(argv[10]);
			#endif
			
		par[i].resultados= crea_array_ent(par[i].NumPix);
		par[i].dix = (int) par[i].anchura/par[i].deltaX;
		par[i].diy= (int) par[i].anchura/par[i].deltaY;

		par[i].angmax=par[i].mlenwid/(2*(par[i].anchura + par[i].d));
		par[i].dist=par[i].anchura/2+par[i].fociii;	
}
	
	double *objeto=malloc(par[0].factordim*par[0].dix*par[0].factordim*par[0].diy*sizeof(*objeto));
	pthread_t threads[NumAdq];
	
	deltax=par[0].factordim*par[0].dix/15;

deltay=par[0].factordim*par[0].diy/15;

radio=deltay/2;


				#ifdef SQUARE
			finp=fopen("matriz_objeto.txt", "w");
		//~ objeto=trataobj(fu1, par[0].factordim,  par[0].dix, par[0].diy, par[0].factordim*par[0].dix/20, par[0].factordim*par[0].diy/20);
					for(a=0; a<(par[0].factordim*par[0].dix); a++)
					{
								for(b=0; b<(par[0].factordim*par[0].diy); b++)
								{		
									*(objeto+a*par[0].factordim*par[0].diy+b)=0.02;
											if(a < par[0].factordim*par[0].dix/2 + deltax && a > par[0].factordim*par[0].dix/2 -deltax && 
											b < par[0].factordim*par[0].diy/2 +deltay && b > par[0].factordim*par[0].dix/2 - deltay) /*Aqui parentesis en las condiciones ¿? */
												*(objeto+a*par[0].factordim*par[0].diy+b)=5;
														//~ printf("Flor");
														//~ getchar();}
									fprintf(finp, "%lf  ", *(objeto+a*par[0].factordim*par[0].diy+b));
							    }
					    fprintf(finp, "\n");
					}
			fclose(finp);
		//~ printf("Dimensiones X = %d, dimension y = %d, deltax=%lf; deltay =%lf\n", par[0].factordim*par[0].dix, par[0].factordim*par[0].diy, deltax, deltay);
		//~ getchar();
			#endif
	
	#ifdef CIRCLE
	
	
		for(a=0; a<par[0].factordim*par[0].dix; a++)
		{
				
				for(b=0; b<par[0].factordim*par[0].diy; b++)
				{				
					*(objeto + a*par[0].factordim*par[0].diy + b)=0;					
							if(((a-par[0].factordim*par[0].dix/2)*(a-par[0].factordim*par[0].dix/2) + (b - par[0].factordim*par[0].diy/2)*(b - par[0].factordim*par[0].diy/2)) < (radio*radio))
							*(objeto + a*par[0].factordim*par[0].diy + b)=50;
				
				}
			//~ guarda_array(objeto, fac*dix, fac*diy, fui); 		
		//~ fclose(fui);		
		}
		//~ objeto=trataobj2(fu1, par[0].factordim, par[0].dix, par[0].diy, par[0].factordim*par[0].dix/10);
	#endif
	dimtotal=par[0].factordim*par[0].dix*par[0].factordim*par[0].diy;
	#ifdef PHANTOM
	
	for(a=0; a<dimtotal; a++) /*inicializa todo el objeto*/
	{
			*(objeto+a)=0;
	}
	for(a=0; a<NElip; a++)
	{
		elipse(par[0].dix, par[0].diy, par[0].factordim, objeto, val[a], x0[a], y0[a], ax[a], bx[a]);
		
	}
	#endif
	
	#ifdef CIRCLE_IN_SQUARE
			finp=fopen("matriz_objeto.txt", "w");
		for(a=0; a<(par[0].factordim*par[0].dix); a++)
			{
					for(b=0; b<(par[0].factordim*par[0].diy); b++)
						{		
							*(objeto+a*par[0].factordim*par[0].diy+b)=0.02;
							
								if(a < par[0].factordim*par[0].dix/2 + deltax && a > par[0].factordim*par[0].dix/2 -deltax && b < par[0].factordim*par[0].diy/2 +deltay && b > par[0].factordim*par[0].dix/2 - deltay) /*El objeto aquí está mal definido */
								*(objeto+a*par[0].factordim*par[0].diy+b)=1;
								
								if(((a-par[0].factordim*par[0].dix/2)* (a-par[0].factordim*par[0].dix/2) + (b - par[0].factordim*par[0].diy/2)*(b - par[0].factordim*par[0].diy/2)) < (radio*radio))
								*(objeto + a*par[0].factordim*par[0].diy + b)=5;
								
								fprintf(finp, "%lf  ", *(objeto+a*par[0].factordim*par[0].diy+b));
					    }
					    fprintf(finp, "\n");
			}
	fclose(finp);
	
	#endif
	
comienzo=clock();



for(i=0; i<NumAdq; i++)
{
	
	par[i].obj=objeto;
	
}
//~ printf("Numero maximo de threads =%d\n",  omp_get_max_threads());
//~ getchar();
omp_set_num_threads(NumAdq);

printf("Numero maximo de threads =%d y NumAng =%d\n",  omp_get_num_threads(), NumAdq);
getchar();
#pragma omp parallel  
{
	#pragma omp for private(cont)
			for(cont=0; cont < NumAdq; cont++) /*Iteracion de angulos */
			{
			
			//~ printf("Numero de hilos = %d\n", omp_get_max_threads());
			
				par[cont].contador=cont;
				#ifdef DOSPI
				par[cont].beta=deltaang*cont;
				#endif
				
				#ifdef ABANICO
				par[cont].beta=angulos[cont];
				printf("Angulo =%lf\n", angulos[cont]);
				#endif
					//~ beta=par[cont].beta;
			
			nth+=1;
			propagacion(&par[cont]);
			
			
				 /*vuelvo a inicializar el array resultados*/
			
				 
				 //~ fu2=fopen("Cuenta de tiempo.txt", "a+");
				
			printf("Final del ángulo %lf\n", par[cont].beta);
				 
			}
			}
	 
	 //~ for (i=0; i<NumAdq; ++i) /*Este for dónde va ¿?¿?*/
	 //~ {
    //~ rc = pthread_join(threads[i], NULL);
		//~ }

fin=clock();
printf("El tiempo de ejecución es %lf y el numero de threads es %d\n", (double)(fin-comienzo)/CLOCKS_PER_SEC, nth);


	return 0;
}



double *crea_array(int di)
{

	int a;
    double *p =malloc(di*sizeof *p);

	for(a=0; a<di; a++)
	{
	*(p+a)=0;
}
	return p;



}

 int *crea_array_ent(int dilm)
{

	int a;
 int *p = malloc(dilm*sizeof (*p)); /*Error de memoria aquí*/
//int p[dilm];

	for(a=0; a<dilm; a++)
	*(p+a)=0;

	return p;

}
void trataobj(FILE *fui, int fac, int dix, int diy, int deltax, int deltay)
{
	
int alm =fac*fac*dix*diy;
int a,b,facdix, facdiy;
double objeto[alm];
facdix=fac*dix;
facdiy=fac*diy;
fui=fopen("Matriz_objeto2.txt", "w+");

	

	for(a=0; a<facdix; a++)
	{
		for(b=0; b<facdiy; b++)
		{		
			*(objeto+a*facdiy+b)=0.01;
		if(a < facdix/2 + deltax && a > facdix/2 -deltax && b < facdiy/2 +deltay && b > facdiy/2 - deltay) /*El objeto aquí está mal definido */
		*(objeto+a*facdiy+b)=5000;
	    }
	}

	guarda_array(objeto, fac*dix, fac*diy, fui); /*Una vez guardado aqui se podría leer desde aquí*/

fclose(fui);
	 
	 


}
double *prop(double delta, double *po1, double *po0) /*Propagacion por el material*/

{





		*po1=*po0+*(po0+1)*delta; /*No creo que sea tan fácil*/
		*(po1+1)=*(po0+1);
		*po0=*po1;
		*(po0+1)=*(po1+1);

		return po0;





}
void guarda_array(double *array, int dimx, int dimy, FILE *fau)
{
	int i, j;

double r;
	for(i=0; i<dimx; i++)

{
for(j=0; j<dimy; j++)
{

 r=*(array + dimx*i + j);
 fprintf(fau, "%lf ", r);


}
fprintf(fau, "\n"); /*Nueva modificación*/

}

/*fclose(fau);*/
}

void guarda_resultados(int *arrayn, int npi,  FILE* finp)
{
int cont;
int res;



for(cont=0; cont<npi; cont++)
{
res=*(arrayn + cont);
fprintf(finp, "%d\n",  res);


}

}

void guarda_resul(double** array, int NumRayos, FILE  *finp, double dist1, double dist2)
{
	
	int cont, cont2;
	double corx;
	corx = dist1;

	for(cont=0; cont<2; cont++)
	{
		for(cont2=0; cont2<NumRayos; cont2++)
		{
			fprintf(finp, "%lf %lf", corx, array[cont2][cont]);
				fprintf(finp, "\n");
		}
	
		corx=dist2;
	}
	
	
}

int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet)
{
//~ int *pix=malloc(2*sizeof(int));
//~ int *pix2=malloc(2*sizeof(int));

static int pix[2];
int pix2[2];

		*pix=rint(dix/2-j*deltaPas/deltX-1);
		*(pix+1)=rint(((*pos-anchura_mall*factodim/2)/deltY));
		*pix2=*pix*cos(bet)+*(pix+1)*sin(bet);
		*(pix2+1)=*(pix+1)*cos(bet)-*(pix)*sin(bet);
		*pix=*pix2+factodim/2*dix;
		*(pix+1)=*(pix2+1)+factodim/2*diy;

return pix;
}


FILE *apertura_archivo(void)
{

char nombrearchivo[25]="FPwML.txt";
//~ #ifdef DOSPI
//~ sprintf(nombrearchivo, "angulo_%d.txt", contador);
//~ #endif

//~ #ifdef ABANICO
//~ sprintf(nombrearchivo, "angulo_bundle%d.txt", contador);
//~ #endif
//~ printf("Se ha abierto el archivo %d\n", contador);

//~ printf("Direccion de memoria del archivo %d nomás abierto=%p", contador, fopen(nombrearchivo,"w"));


//~ if(fopen(nombrearchivo, "w") == 0x000000);
//~ {
	//~ printf("Error en %d", contador);
	

	
//~ }
return fopen(nombrearchivo, "w");



}

double *calcula_mu (int longitud, double *mul, FILE *fin) /*Esta subrutina añade una centesima de segundo por ejecucion. En total 0.01*NumAng*Niter*N*deltaPaso/deltaX segundos (mas o menos la edad del Universo) */
{
int contador;
double res;
clock_t start, stop;
start=clock();
for(contador=0; contador < longitud;  contador++)
{

fscanf(fin, "%lf", &res);

}
rewind(fin);
*mul=res;
stop=clock();
/*printf("La subrutina tarda %lf\n", ((double) (stop-start)/CLOCKS_PER_SEC));*/

return mul;

}

void *propagacion(void *threadarg) /*Falta lo del objeto --- Pasar por el struct los parametros de la main lens y de las lenslet*/
{
	double num_rand, betarnd, deltadis, deltaP, Niter,  mu, cosang, sinang, f , focii, mlenswidth, distancia;
	double pos0[2], pos1[2];
	int i, j, pixx, pixy, num_pix, NRayos, cont, cont2, deltaang, pi1, pi2, pami, pixu, numlens, mlensdiv, ind;
	int *pix=malloc(2*sizeof(int));
	int *pix2=malloc(2*sizeof(int));
	struct param *mispar;
	mispar=(struct param *) threadarg;
	cont = mispar -> contador;
	numlens = mispar -> nulen;
	mlensdiv= mispar -> mlendiv;
	mlenswidth= mispar -> mlenwid;
	f = mispar -> ff;
	focii = mispar -> fociii;
	distancia = mispar -> dist;
	double *objet = mispar -> obj;
	int NumPi = mispar->NumPix;
	double anchura = mispar->anchura;
	double angma = mispar -> angmax;
	double deltaPas = mispar ->deltaPaso;
	int dix =mispar-> dix;
	int diy=mispar ->diy;
	double deltX = mispar -> deltaX;
	double deltY = mispar -> deltaY;
	int factor = mispar -> factordim;
	double bet = mispar -> beta;
	 int *resultados=malloc(NumPi*sizeof(int)); /*Error de memoria aquí*/
	resultados=crea_array_ent(NumPi);
	
	//~ printf("Delta Y = %lf, deltaX = %lf, beta = %lf, dimx=%d", deltX, deltY, bet, dix);

	//~ int *resultados = malloc(NumPi*sizeof(unsigned long int));
	//~ resultados=crea_array_ent(NumPi); /*La creación aquí del array resultados petará posteriormente*/
	//~ deltaPaso=0.001;
	FILE *fui, *fout;
	fui=fopen("Matriz_objeto2.txt", "r+");
	NRayos = mispar -> N;
	deltaP=anchura/NumPi; /*Tamaño de cada detector*/

srand(time(NULL));

Niter=(int)(dix*deltX/deltaPas);
	double deltadi = mispar -> deltadis;
	double delta = mispar -> d;	
	fout=fopen("Rayos.txt", "w+");
	printf("Contador = %d\n", cont);
	cosang=cos(bet);
	sinang=sin(bet);
	pami=0;
	double resul[NRayos][2]; /*Inicializarlo a -1 ¿?*/

			for(cont=0; cont<2; cont++)
			{
				for(cont2=0; cont2<NRayos; cont2++)
				{
				resul[cont2][cont]=-1;
				}
			}
			ind=0;
			
for(i=0; i<NRayos; i++) /*Iteracion de todos los rayos para un mismo angulo*/
{

#ifdef RAND_ANG
	num_rand=rand()/((double)RAND_MAX +1);
	betarnd=-angma+2*angma*num_rand;
	num_rand=rand()/((double)RAND_MAX +1);
	*pos0 = factor*anchura/2-deltadis+ 2*deltadis*num_rand; /*Fuente extensa, posicion y aleatoria*/
	 *(pos0+1)= betarnd;
	 
	 #endif
	 
	 #ifdef DEB_ANG /*angulos equiespaciados en el cono de posibles angs*/
	 num_rand=2*angma/NRayos*i; /*Este factor 2 ¿?*/
	 betarnd=-angma+2*angma*num_rand;
	num_rand=0.5;
	*pos0 = factor*anchura/2-deltadis+ 2*deltadis*num_rand; /*Fuente extensa, posicion y aleatoria*/
	 *(pos0+1)= betarnd;
	 #endif
	 
	*pos1=*(pos0+1)*delta+*pos0;
	*pos0=*pos1;
	*(pos1+1)=*(pos0 +1 );

fprintf(fout, "%f %f %f %d %d\n", pos0[0], pos0[1],betarnd, i, 0);
	pixx=-1;
	pixy=-1;



	for(j=0; j<Niter; j++) /*Recorrido de cada rayo*/

		{

		*(pos1)=*(pos0+1)*deltaPas+*pos0;
		*(pos0)=*(pos1);
										//~ pix=my_pixel(dix, diy, j, deltaPas, deltX, deltY, factor, anchura, pos0, bet); /*aquí se calcula un numero de pix coherente para distintos ángulos de adquisicion*/

		*pix=rint(dix/2-j*deltaPas/deltX-1);
		*(pix+1)=rint(((*pos0-anchura*factor/2)/deltY));
		*pix2=*pix*cosang+*(pix+1)*sinang;
		*(pix2+1)=*(pix+1)*cosang-*(pix)*sinang;
		*pix=*pix2+factor/2*dix;
		*(pix+1)=*(pix2+1)+factor/2*diy;
	
								    /*fseek(fu1, (pixx*dimy*factordim+pixy)*sizeof(double), SEEK_SET);  he quitado esto (pixx*dimy*factordim+pixy)*
									fscanf(fu1, "%lf", &mu);	Atenuacion lineal
									rewind(fu1);*/
	
		pi1=*pix;
		pi2=*(pix +1);
		mu=*(objet + pi1*diy*factor+pi2);
		
					if(pixx != *(pix) || pixy != *(pix+1)) /*Esto es un OR */
								{
									pixx=*(pix);
									pixy=*(pix+1);
									//~ calcula_mu((pixx*diy*factor+pixy), &mu, fui); /*Very time-consuming*/
									num_rand=rand()/((double)RAND_MAX +1);
											
													if(num_rand > exp(-mu*deltX)) 
															{//~ printf("Rayo no llego en angulo %d para mu =%lf\n", cont, mu);
																pami+=1;
																break;
															}
					
								}
				
					//~ if(j== (Niter-1))
						//~ {
						    //~ num_pix=rint((*pos0-factor*anchura/2+anchura/2)/deltaP); /*Num de detector al que llega el rayo*/							
							//~ printf("El número de píxel es %d y pos del rayo es %lf, mu =%lf y angulo %d, mientras, no han llegado %d rayos\n", num_pix, *pos0-factor*anchura/2+anchura/2, mu, cont, pami);
							//~ pami=0;
													//~ if(num_pix >= 0 && num_pix < NumPi)
															//~ *(resultados+num_pix)=*(resultados+num_pix)+1;
					
						//~ }
				

} /*Final de la interacción de los X-Ray a través del objeto*/

/*Poner aquí el if(j==(Niter-!)) --- probablemente*/
/* Propagation through main lens  */
/*Realizar una traslación para la pos0*/
if(j==(Niter-1))
{printf("Llego al if y pami = %d, numero de rayo %d\n", pami, i);
*pos0=*pos0 + (1-factor)*anchura/2; /*Traslacion--- revisar*/

*pos1=*pos0+*(pos0+1)*distancia;
*(pos1+1)=-*(pos0 +1);
*pos0=*pos1;
*(pos0+1)=*(pos1+1);



pixu=rint(*pos0/mlenswidth*mlensdiv);/* Calcular aquí en qué punto de la lente corta*/
*pos1=*pos0+(*pos0+1)*focii; /*Propagación hasta la minilente*/
*pos0=*pos1;
resul[ind][0]=pixu*mlenswidth/mlensdiv;
pixu=rint(*pos0/mlenswidth*numlens);/*Calcular aqui la lente a la que llega */

*pos0=pixu*mlenswidth/numlens;/*igualar la pos0 a la posicion central de la lentilla*/

*pos1=*pos0 + *(pos0+1)*f;

pixu=rint(*pos1/deltaP);/*Calcular pixel al que llega con pos1*/
resul[ind][1]=pixu*deltaP; /* Guardarlo no como num de pixel sino como  coordenada y*/
	pami=0;
	ind+=1;
}
}
fclose(fui);
fclose(fout);

	 fui=apertura_archivo(); /*Crea un archivo para cada ángulo*/

	guarda_resul(resul, NRayos, fui, (distancia + anchura),  (distancia + anchura + focii + f)); /*Guarda los rayos detectados por cada pixel*/
	fclose(fui);
	free(resultados);
	

return NULL;

}
