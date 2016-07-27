#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
/*Inicializar las funciones*/
double *prop(double delta, double *po1, double *po0);
double *crea_array(int dd);
void guarda_array(double *arra, int dix, int diy, FILE *fuuu);
double *crea_matriz(int dx, int dy, int dex, int dey);
unsigned long int *crea_array_ent(int dil);
void guarda_resultados(unsigned long int *arra, int nmpix,  FILE *fiiii);
FILE *apertura_archivo(int con);
int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet);
/*double *matrixmul(double aaa, double bbb);*/
/*#define anchura_malla 1
#define deltaY 0.001
#define deltaX 0.001
#define dimx anchura_malla/deltaX
#define dimy anchura_malla/deltaY*/
#define PI 3.141592653
int main(int argc, char **argv)
{
int i, j,  cont,  hit, no_hit, dimx, dimy, num_pix, factordim, NumAdq, Niter, pixx, pixy;
double  d, deltaX, deltaY, deltaP, betarnd, num_rand, mu, anchura_malla, delta, deltadis, deltaang, beta, deltaPaso;
unsigned long int N, NumPix;
factordim=2; /*Numero de veces que es más grande el grid mesh utilizado -siempre cuadrado- que el objeto a estudiar -también cuadrado-*/
NumAdq=20; /*Numero de angulos de adquisición*/
deltaang=2*PI/NumAdq;
anchura_malla=1; /*En unidades de medida, anchura del objeto*/
deltadis=0.001;
FILE  *fu1, *fprueba;
d=0.5;
N=100000; /*Número de rayos por cada ángulo de adquisición*/
deltaX=deltaY=0.01; /*Intervalo o paso del mesh grid*/
deltaP=0.01; /*Tamaño de cada detector*/
deltaPaso=0.001; /*Intervalo de propagación del rayo en cada iteración*/
delta=d; /*Distancia de la fuente al objeto*/
hit=0;
no_hit=0;
NumPix=(int)anchura_malla/deltaP; /*Numero de detectores*/
dimx=(int)(anchura_malla/deltaX); /*Num de nodos en direccion x*/
dimy=(int)(anchura_malla/deltaY);
Niter=(int)(dimx*deltaX/deltaPaso);/*Numero de iteraciones/propagaciones de cada rayo*/


double *pos0 = malloc(2*sizeof(double));

double *pos1 = malloc(2*sizeof(double));

int *pix =malloc(2*sizeof(int));
double *objeto = malloc(factordim*factordim*dimx*dimy*sizeof(double));
unsigned long int *resultados = malloc(NumPix*sizeof(unsigned long int)); /*Igual esto no necesario*/
objeto=crea_matriz(factordim*dimx, factordim*dimy, factordim*dimx/10, factordim*dimy/10);


 /*Otra funcion crea_array con retorno de entero*/
double angmax=anchura_malla/(2*(anchura_malla+d)); /* [-angmax, angmax] es el intervalo angular en el que proyecto rayos*/

srand(time(NULL));

printf("Comienzo proyeccion con varios angulos\n");


for(i=0; i<factordim*dimx; i++)
{for(j=0; j<factordim*dimy; j++)
	{
if(*(objeto+i*dimx+j) != 0)
	printf("%lf\n",*(objeto +i*dimx*factordim+j));
	printf("Da el fallo en el elemento %d, %d; el elemento es: %lf \n", j, i,*(objeto +i*dimx+j));
	}

}

getchar();

fu1=fopen("Matriz_objeto.txt", "w+");
guarda_array(objeto, factordim*dimx, factordim*dimy, fu1); /*Una vez guardado aqui se podría leer desde aquí*/
free(objeto);
fclose(fu1);
fu1=fopen("Matriz_objeto.txt", "r+");



/*freopen() ¿?¿?¿?*/
for(cont=0; cont < NumAdq; ++cont) /*Iteracion de angulos */
{
beta=cont*deltaang;
for(i=0; i<N; ++i) /*Iteracion de todos los rayos para un mismo angulo*/
{
	num_rand=rand()/((double)RAND_MAX +1);
	betarnd=-angmax+2*angmax*num_rand;
num_rand=rand()/((double)RAND_MAX +1);
	*pos0 = factordim*anchura_malla/2-deltadis+ 2*deltadis*num_rand; /*Fuente extensa, posicion y aleatoria*/
	 *(pos0+1)= betarnd;



	*pos1=*(pos0+1)*delta+*pos0;
	*pos0=*pos1;
	*(pos1+1)=*(pos0 +1 );
/*fprintf(fout, "%f %f %f %d %d\n", pos0[0], pos0[1],betarnd, i, 0);*/

resultados=crea_array_ent(NumPix); /*Aquí debe de estar el problema de memory allocation*/
	for(j=0; j<Niter;++j) /*Recorrido de cada rayo*/

	{

	*pos1=*(pos0+1)*deltaPaso+*pos0;
	*pos0=*pos1;
/*Aquí va my_pixel */
pix=my_pixel(dimx, dimy, j, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0, beta); /*aquí se calcula un numero de pix coherente para distintos ángulos de adquisicion*/
pixx=*pix;
pixy=*(pix+1);
/*  mu sacarlo del archivo */

    fseek(fu1, (pixx*dimy*factordim+pixy)*sizeof(double), SEEK_SET); /* he quitado esto (pixx*dimy*factordim+pixy)* */
	fscanf(fu1, "%lf", &mu);	/* Atenuacion lineal*/
	rewind(fu1);
printf("La mu de %d %d es %lf\n",pixx, pixy, mu); /*Da muchos ceros */

		num_rand=rand()/((double)RAND_MAX +1);

      /*  printf("el angulo es %d la N es = %dla j es = %d y el pixel %d %d\n", cont, i,  j, *pix, *(pix+1));*/




		if(num_rand > exp(-mu*deltaX)) /*Pasa o no el rayo*/
		{no_hit=no_hit+1;
printf("Rayo no llego\n");
			break;
	}

}
if(j== Niter)
{
	num_pix=rint((*pos0-factordim*anchura_malla/2+anchura_malla/2)/deltaP); /*Num de detector al que llega el rayo*/
/*printf("numpix = %d, pos0 = %lf, pos0prim = %lf beta = %lf\n", num_pix, *pos0, *pos0-factordim*anchura_malla/2+anchura_malla/2, beta);*/

if(num_pix >= 0 && num_pix < NumPix)
	*(resultados+num_pix)=*(resultados+num_pix)+1;
printf("Rayo llego\n");
	hit=hit+1;
	}
}
printf("Final del ángulo %lf", beta);
fprueba=apertura_archivo(cont); /*Crea un archivo para cada ángulo*/

	guarda_resultados(resultados, NumPix, fprueba); /*Guarda los rayos detectados por cada pixel*/
	fclose(fprueba);
	free(resultados);
	 /*vuelvo a inicializar el array resultados*/
}
printf( "\nHan llegado = %d \nNo han llegado %d de %lu rayos\n", hit, no_hit, N*NumAdq);





	return 0;
}


double *crea_matriz(int dix, int diy, int deltax, int deltay) /*Esto no se libera */
{
int a, b;
double *r=malloc(dix*diy*sizeof(*r));
	for(a=0; a<dix; a++)

	{
		for(b=0; b<diy; b++)

		{		*(r+a*diy+b)=10;

		if(a < dix/2 + deltax && a > dix/2 -deltax && b < diy/2 +deltay && b > diy/2 - deltay) /*El objeto aquí está mal definido */

		*(r+a*diy+b)=10 ;






	}

		}

		return r;



}
double *crea_array(int di)
{

	int a;
    /*double *p =malloc(di*sizeof *p);*/
double *p;
	for(a=0; a<di; a++)
	{
	*(p+a)=0;
}
	return p;



}

unsigned long int *crea_array_ent(int dilm)
{

	int a;
	  unsigned long int *p = malloc(dilm*sizeof (*p));

	for(a=0; a<dilm; a++)
	*(p+a)=0;
free(p);
	return p;

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


}

/*fclose(fau);*/
}

void guarda_resultados(unsigned long int *arrayn, int npi,  FILE* finp)
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
int *pix=malloc(2*sizeof(int));
int *pix2=malloc(2*sizeof(int));

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

return fopen(nombrearchivo, "w");



}
