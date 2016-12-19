#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "funciones.h"
#include <omp.h>
#include <errno.h>
#define deltaX 0.02
#define deltaY 0.02
//~ #define NumPix 100 /*Calcular numpix y numang automaticamente ¿?¿?*/
/*#define RANDINPUT*/
//~ #define anchura_malla 1
#define INPUT_TXT
#define INPUT_ARCHIVO
#define PI 3.141592653
#define RAND_ANG
#define SUM_ONE /*Suma un entero al pixel correspondiente en la matriz resultante*/
#define DOSPI
void guarda_matriz(int *res, int dix, int diy, FILE *fghq);
double *import_data(double *datain);
int *crea_matriz(int dim);
unsigned long int *crea_array (int ded);
double *my_pixel2(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet);
//~ void *crea_imagen( int ray, int numan);
void my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double sinbet, double cosbet, int *numpix);
FILE *apertura_archivo(int co);
FILE *apertura_archivo2(int co, int co2, int co3);
int main(int argc, char **argv)
{

int dimx, dimy, ray, i, cont,  j, q, pasx, factordim, NumAng, numray, conta, Num_Pix, dimtotal; /*Variables n, d no usadas*/

int array_inp[4];
double  yp, yl, division,  f, betarand, betamin, betamax, ang, beta, deltaPaso, deltaL, deltaP, anchura_malla, tiempo, cosbeta, sinbeta;
double pos0rnd[2], pos1rnd[2], angulos[5];
clock_t comienzo, fin;
FILE *finp, *fout, *fpru;

 //~ int *numpix=malloc(2*sizeof(int));
 //~ #ifdef SUM_ONE
//~ int *numpix=malloc(2*sizeof(int)); /*No bastaria con numpix[2]*/
//~ #endif

//~ #ifdef WEIGHT_SUM

//~ double *numpix = malloc(4*sizeof(double));

//~ #endif
//~ static int numpix[2];



 #ifdef INPUT_TXT
finp=fopen("inputint.txt", "r");

for(cont=0; cont<4; cont++)
{
fscanf(finp, "%d", &array_inp[cont]);
printf("array_input[%d]=%d\n",cont, array_inp[cont]);


}
fclose(finp);
factordim=array_inp[3];
Num_Pix=array_inp[2];
NumAng=array_inp[0];

#endif 

#ifdef ABANICO 

finp=fopen("NumAngAbanico.txt", "r");
fscanf(finp, "%d", &NumAng);
fclose(finp);
finp=fopen("AngulosBundle.txt", "r");
for(cont=0; cont<NumAng; cont++)
{
	fscanf(finp, "%lf", &angulos[cont]);
	printf("%lf\n", angulos[cont]);
	
}
fclose(finp);
//~ getchar();
#endif

anchura_malla=1;

conta=0;

ray=0;
ang=2*PI/NumAng;
dimx=(int)anchura_malla/deltaX;
dimy=(int)anchura_malla/deltaY;
deltaP=anchura_malla/Num_Pix;
deltaL=deltaP*5;
deltaPaso=0.0001; /*Intervalo de propagacion de cada rayo*/
 /*Factor dim se deberia meter de output*/
pasx=(int)dimx*deltaX/deltaPaso;

dimtotal=factordim*factordim*dimx*dimy;
//~ int *matriz_resul = malloc(dimtotal*sizeof(*matriz_resul));/*Esta función reserva 10 Mb de memoria ---- utilizar la funcion calloc(dimtotal, sizeof(int)) ¿?¿?¿*/
//~ unsigned long int *datainp = malloc(Num_Pix*sizeof(*datainp));  /*Reserva de memoria */ 
//~ matriz_resul=crea_matriz(dimtotal); /* Crea el mallado que va a guardar los resultados*/
//~ datainp=crea_array(Num_Pix);
#ifdef SUM_ONE
//~ int *matriz_resul=calloc(dimtotal, sizeof(*matriz_resul));

int (*matriz_resul)[factordim*dimy]=malloc(sizeof(int[factordim*dimx][factordim*dimy]));


for(cont=0; cont<(factordim*dimx); cont++)
{
	for(i=0; i<(factordim*dimy); i++)
	{
		
		matriz_resul[cont][i]=0;
	}
}
#endif


#ifdef WEIGHT_SUM
double *matriz_resul=calloc(dimtotal, sizeof(*matriz_resul));
#endif

unsigned long int (*datainp)[NumAng] = malloc(sizeof(unsigned long int[Num_Pix][NumAng]));
//~ omp_set_num_threads(NumAng);
#ifdef INPUT_ARCHIVO

	for(cont=0; cont<NumAng; cont++)
	{finp=apertura_archivo(cont); /*AQUI IFDEF PARA DOSPI/ABANICO*/
		for(i=0; i<Num_Pix; i++)
		{
		fscanf(finp, "%lu", &datainp[i][cont]);
		}
	fclose(finp);
	}

#endif

#ifdef RANDINPUT /*Inicializacion de un array con num enteros aleatorios*/
	srand(time(NULL));
	for(cont=0; cont<NumAng; cont++)
	{
		for(i=0; i<Num_Pix; i++)
		{
		datainp[i][cont]=rand() % 11;
		
		}
	}
#endif
comienzo=clock();
 //~ #ifdef SUM_ONE
//~ int *numpix=malloc(2*sizeof(int)); /*No bastaria con numpix[2]*/
//~ static int *numpix;
//~ #endif

#ifdef WEIGHT_SUM

double *numpix = malloc(4*sizeof(double));

#endif
//~ omp_set_num_threads(NumAng);
#pragma omp parallel 
{
	#ifdef SUM_ONE
		int numpix[2];
	#endif
	
	//~ #pragma omp for shared(datainp, matriz_resul) private(j, cont,q, beta, i, numray,sinbeta, cosbeta) 
	
	#pragma omp for private(j, cont,q, beta, i, numray,sinbeta, cosbeta) 
	for(cont=0; cont<NumAng; cont++) /*Cambiar este for de lugar ya que solo esta dentro del INPUT_ARCHIVO ¿?*/
	{
		#ifdef DOSPI
		beta=cont*ang;
		#endif
		#ifdef ABANICO
		beta=angulos[cont];
		#endif
		cosbeta=cos(beta);
		sinbeta=sin(beta);
		int numpixx2, numpixy2;
		
		
		printf("Comienza el cuerpo de la simulación, beta=%lf, cosbeta=%lf, sinbeta=%lf, cont=%d, num de hilo=%d, NumAng=%d\n", beta, cosbeta, sinbeta, cont, omp_get_thread_num(), NumAng);
		for(i=0; i<Num_Pix; i++)/*Cuerpo de toda la simulación, proyeccion de cada pixel */
		{
			//~ printf("numpix=%d, hilo=%d\n", i, omp_get_thread_num());
			if(datainp[i][cont] != 0)
			{
				/*sabiendo la posicion del detector y la distancia lenslet-detector se tiene un cono angular dentro del cual se escoge un angulo aleatoriamente*/
				yp=i*deltaP+deltaP/2;
				division=yp/deltaL;
				division=(int)division;
				yl=division*deltaL+deltaL/2;
				//~ d=f=0.0001;
				
				f=0.7;
				betamin=(yp-yl-0.5*deltaP)/f;
				betamax=(yp-yl+0.5*deltaP)/f;
				
				numray=datainp[i][cont];
				for(j=0; j<numray; j++)
				{
					ray=ray+1;
					#ifdef RAND_ANG
						betarand=rand()/((double) RAND_MAX +1)*(betamax-betamin)+betamin; 
					#endif
					#ifdef ANG_DEB
						betarand=0;
					#endif
					/*betarand=0.05;*/
					pos0rnd[0]=yl;
					pos0rnd[1]=betarand;
					numpixx2=-1;
					numpixy2=-1; /* Propagacion de delta ¿?  */
					//~ fpru=apertura_archivo2(i, j, cont);
					//~ fprintf(fpru,  "%lf   %lf	%d		%d\n", pos0rnd[0], pos0rnd[1], numpixx2, numpixy2);
					for(q=0; q<pasx; q++)
					{
						pos1rnd[0]=pos0rnd[0]+pos0rnd[1]*deltaPaso;
						pos1rnd[1]=pos0rnd[1];
						pos0rnd[0]=pos1rnd[0];
						
						#ifdef SUM_ONE
						my_pixel(dimx, dimy, q, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0rnd, sinbeta, cosbeta, numpix); 
						//~ printf("numpix=%d, numpix2=%d\n", *numpix, *(numpix+1));
						#endif
						#ifdef WEIGHT_SUM
						numpix=my_pixel2(dimx, dimy, q, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0rnd, beta); 
						#endif
						
						if(numpixx2!=*numpix || numpixy2!=*(numpix+1)) /*Este if probablemente sobre pero me libra de posibles abortos del programa por error */
						{
							numpixx2=*numpix;
							numpixy2=*(numpix+1);
							//~ numpixx2=2;
							//~ numpixy2=4;
							//~ printf("numpix2 = %d, numpiy2=%d, numpix=%d, numpix+1=%d, numdepixel=%d,numderayo=%d, datainp=%d, hilo=%d, cont=%d\n", numpixx2, numpixy2, *numpix, *(numpix+1),i,j, datainp[i][cont], omp_get_thread_num(), cont);
							matriz_resul[numpixx2][numpixy2]+=1;
							conta+=1;
							//~ fprintf(fpru, "%lf   %lf	%d		%d\n", pos0rnd[0], pos1rnd[1], numpixx2, numpixy2);
						}
					}
				}
				//~ fclose(fpru);
			}			
		}
		
		fin=clock();
		//~ tiempo=difftime(fin, comienzo);
		//~ fpru=fopen("Calculo_tiempo.txt", "a+");
		
		//~ fprintf(fpru, "%lf\n", tiempo);
		//~ fclose(fpru);
		
		/*Guarda todos los resultados en una matriz*/
		
	}
	printf("Fin del angulo %d\n", cont);
}

#pragma omp barrier
//~ free(numpix);

fin=clock();
printf("Tiempo de ejecucion de los %d hilos = %lf\n", NumAng, (double)(fin-comienzo)/CLOCKS_PER_SEC );
char archivo[50]="Resultados_BP.txt";
fout=fopen(archivo, "w+");


guarda_matriz(matriz_resul,factordim*dimx, factordim*dimy, fout);
free(matriz_resul);
free(datainp);

fclose(fout);

crea_imagen(ray, NumAng);
//~ FILE *pipeplot=popen("gnuplot -persist", "w");
//~ /*fprintf(pipeplot, "set title /lo que sea/"   */
//~ fprintf(pipeplot, "set pm3d map\n");
//~ fprintf(pipeplot, "set size square 1, 1\n");
//~ fprintf(pipeplot, "splot '/home/pablo/Downloads/Resultados_BP.txt' matrix\n");
//~ sleep(5);
//~ fprintf(pipeplot, "set term png\n");
//~ fprintf(pipeplot, "set output 'BP_Imagen_Angulos%d_Rayos%d.png'\n", NumAng, ray);
//~ fprintf(pipeplot, "replot\n");

//~ fflush(pipeplot);
//~ fprintf(pipeplot, "quit\n");
//~ pclose(pipeplot);
//~ printf("conta=%d\n", conta);
//~ getchar();
return 0;
}





int *crea_matriz(int dim)
{
int a;


int *r =malloc(dim*sizeof(*r)); /*Reserva de memoria unnecessary*/

for(a=0; a<dim; a++)
{
	
	*(r+a)=0;
}
	//~ for(a=0; a<dix; a++)

	//~ 
		//~ for(b=0; b<diy; b++)

		//~ {

		//~ r[a*diy+b]=0;


	//~ }

		//~ free(r);
		return r;


}



double *import_data(double *datai) 
{
size_t sizedata;
sizedata=sizeof(datai)/sizeof(datai[0]);
int i;
double *p = malloc(sizedata*sizeof(*p));
for(i=0; i<sizedata; i++)
{
p[i]=datai[i];


}

return p;

}

unsigned long int *crea_array(int di)
{

	int a;
	  unsigned long int *p = malloc(di*sizeof(*p)); /*Reserva de memoria incorrecta */

	for(a=0; a<di; a++)
	p[a]=0;

	return p;

}

void guarda_matriz(int *resu, int dimxx, int dimyy, FILE *fgu)
{
int con1, con2;

for(con1=0; con1 <dimxx; con1++)
{
for(con2=0; con2 <dimyy; con2++)
{

fprintf(fgu,"%d   ", *(resu+con1*dimxx+con2));
}
fprintf(fgu, "\n");


}
printf("Final de la subrutina guarda_matriz\n");



}
FILE *apertura_archivo(int contador)
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "angulo_%d.txt", contador);

return fopen(nombrearchivo, "r");



}
FILE *apertura_archivo2(int contador, int contador2, int contador3) 
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "rayo%d_%d_%d.txt", contador3, contador, contador2);

return fopen(nombrearchivo, "w");



}


void my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, 
				double anchura_mall, double *pos, double sinbet, double cosbet, int *numpix)
{
	
 
static int pix[2];
int pix2[2];

		*pix=rint((j*deltaPas-anchura_mall/2)/deltX); 
		*(pix+1)=rint(((*pos-anchura_mall/2)/deltY));
		*pix2=*pix*cosbet-*(pix+1)*sinbet;
		*(pix2+1)=*(pix+1)*cosbet+*(pix)*sinbet;
//~ *pix=*pix2+factodim/2*dix;
//~ *(pix+1)=*(pix2+1)+factodim/2*diy;
//~ *numpix=*pix;
//~ *(numpix+1)=*(pix+1);

*numpix=*pix2+factodim/2*dix;
*(numpix+1)=*(pix2+1)+factodim/2*diy;

//~ printf("pix=%d, pix1=%d, *pix2=%d, *(pix2+1)=%d, j=%d, anchura=%lf, factodim=%d, dix=%d\n", *pix, *(pix+1), *pix2, *(pix2+1), j, anchura_mall, factodim, dix);

//~ return pix;


}


double *my_pixel2(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet)
{


static double pix[4];
double pix2[2];

*pix = (j*deltaPas -anchura_mall/2)/deltX;
*(pix+1) =(*pos - anchura_mall/2)/deltY;
*pix2 = *pix*cos(bet)-*(pix+1)*sin(bet);
*(pix2+1)=*(pix+1)*cos(bet)+*(pix)*sin(bet);

*pix =trunc(*pix2);
*(pix+1) =trunc (*(pix2+1));
*(pix + 2) = *pix2 - *pix;
*(pix+3)=*(pix2+1) - *(pix+1);

return pix;

}



