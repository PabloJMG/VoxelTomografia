

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>
#define NTHREADS 10





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
double anchura;
double angmax;


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
FILE *apertura_archivo(int con);
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
#define DEB_ANG
int main(int argc, char **argv)
{
int i, cont, dimx, dimy, factordim, NumAdq,rc;
int array_inp[3];
double  d, deltaX, deltaY, anchura_malla, deltaang, beta,  angmax;
//~ pthread_t threads[NTHREADS];
//~ int thread_arg[NTHREADS];
clock_t comienzo, fin;
FILE  *fu1;



#ifdef INPUT_TXT
fu1=fopen("input.txt", "r");

for(cont=0; cont<3; cont++)
{
fscanf(fu1, "%d", &array_inp[cont]);
printf("array_input[%d]=%d\n",cont, array_inp[cont]);
}
fclose(fu1);
#endif

factordim=2; 


NumAdq=array_inp[0]; /*Numero de angulos de adquisición*/
//~ N=array_inp[1];
//~ NumPix=array_inp[2];
struct param par[NumAdq];
//~ int *resultados = malloc(NumPix*sizeof(unsigned long int));

deltaang=2*PI/NumAdq;
anchura_malla=1; /*En unidades de medida, anchura del objeto*/

deltaX=deltaY=0.02; /*Intervalo o paso del mesh grid*/
 /*Intervalo de propagación del rayo en cada iteración*//*Distancia de la fuente al objeto*/

/*Tamaño de cada detector*/
dimx=(int)(anchura_malla/deltaX); /*Num de nodos en direccion x*/
dimy=(int)(anchura_malla/deltaY);
d=0.05;

/*Esta malloc...*/
 
srand(time(NULL));


trataobj(fu1, factordim,  dimx, dimy, factordim*dimx/6, factordim*dimy/6); /*Crea objeto y lo guarda en un archivo de texto*/
//~ fu1=fopen("Matriz_objeto2.txt", "r+");

angmax=anchura_malla/(2*(anchura_malla+d)); 

	for(i=0; i<NumAdq; i++){
	#ifdef INPUT_TXT
	par[i].N=array_inp[1];
	par[i].NumPix=array_inp[2];
	#endif
	
	#ifdef ARGS
	par[i].N=atoi(argv[2]);
	par[i].NumPix=atoi(argv[3]);
	
	#endif
	par[i].resultados= crea_array_ent(par[i].NumPix);
	par[i].dix =dimx;
	par[i].diy=dimy;
	par[i].deltaPaso=0.001;
	par[i].deltaX=0.01;
	par[i].deltaY=0.01;
	par[i].factordim=2;
	par[i].anchura=anchura_malla;
	par[i].angmax=angmax;
	
}
	
	pthread_t threads[NumAdq];

	
comienzo=clock();
#pragma omp parallel 
#pragma omp for
for(cont=0; cont < NumAdq; cont++) /*Iteracion de angulos */
{
	





	par[cont].contador=cont;
	par[cont].beta=deltaang*cont;
	beta=par[cont].beta;
	//~ thread_arg[i]=i;
	printf("Contador = %d\n", cont);
	

propagacion(&par[cont]);


	 /*vuelvo a inicializar el array resultados*/

	 
	 //~ fu2=fopen("Cuenta de tiempo.txt", "a+");
	
printf("Final del ángulo %lf\n", beta);
	 
}

	 
	 //~ for (i=0; i<NumAdq; ++i) /*Este for dónde va ¿?¿?*/
	 //~ {
    //~ rc = pthread_join(threads[i], NULL);
		//~ }

fin=clock();
printf("El tiempo de ejecución es %lf\n", (double)(fin-comienzo)/CLOCKS_PER_SEC);


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


FILE *apertura_archivo(int contador)
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "angulo_%d.txt", contador);
printf("Se ha abierto el archivo %d", contador);


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

void *propagacion(void *threadarg) /*N=argv[1]*/
{
	double num_rand, betarnd, deltadis, deltaP, Niter, delta, d, mu;
	double pos0[2], pos1[2];
	int i, j, pixx, pixy, num_pix, NRayos, cont, deltaang;
	int *pix=malloc(2*sizeof(int));
	struct param *mispar;
	mispar=(struct param *) threadarg;
	cont = mispar -> contador;
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
	printf("Delta Y = %lf, deltaX = %lf, beta = %lf, dimx=%d", deltX, deltY, bet, dix);

	//~ int *resultados = malloc(NumPi*sizeof(unsigned long int));
	//~ resultados=crea_array_ent(NumPi); /*La creación aquí del array resultados petará posteriormente*/
	//~ deltaPaso=0.001;
	FILE *fui, *fout;
	fui=fopen("Matriz_objeto2.txt", "r+");
	NRayos = mispar -> N;
	deltaP=anchura/NumPi; /*Tamaño de cada detector*/

srand(time(NULL));

Niter=(int)(dix*deltX/deltaPas);
delta=d=0.5;
deltadis=0.0001;	
fout=fopen("Rayos.txt", "w+");


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
	pix=my_pixel(dix, diy, j, deltaPas, deltX, deltY, factor, anchura, pos0, bet); /*aquí se calcula un numero de pix coherente para distintos ángulos de adquisicion*/


    /*fseek(fu1, (pixx*dimy*factordim+pixy)*sizeof(double), SEEK_SET);  he quitado esto (pixx*dimy*factordim+pixy)*
	fscanf(fu1, "%lf", &mu);	Atenuacion lineal
	rewind(fu1);*/
	
if(pixx != *(pix) || pixy != *(pix+1)) /*Esto es un OR */
			{
	pixx=*(pix);
	pixy=*(pix+1);
	calcula_mu((pixx*diy*factor+pixy), &mu, fui);
	num_rand=rand()/((double)RAND_MAX +1);


		if(num_rand > exp(-mu*deltX)) 
				{	
					printf("Rayo no llego\n");
					//~ getchar();
				
					break;
				}

			}

if(j== (Niter-1))
	{
	num_pix=rint((*pos0-factor*anchura/2+anchura/2)/deltaP); /*Num de detector al que llega el rayo*/

printf("El número de píxel es %d y pos del rayo es %lf\n", num_pix, *pos0-factor*anchura/2+anchura/2);
if(num_pix >= 0 && num_pix < NumPi)
			{
				*(resultados+num_pix)=*(resultados+num_pix)+1;

			}

	


		}


}


	
}
fclose(fui);
fclose(fout);
	 fui=apertura_archivo(cont); /*Crea un archivo para cada ángulo*/

	guarda_resultados(resultados, NumPi, fui); /*Guarda los rayos detectados por cada pixel*/
	fclose(fui);
	free(resultados);
	

return NULL;

}
