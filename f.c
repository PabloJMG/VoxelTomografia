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
void trataobj(FILE *fui, int fac, int fac2, int dix, int diy, int deltax, int deltay); /*El fac2 es solo para este tip de objeto*/
//~ double *crea_matriz(int dx, int dy, int dex, int dey);
 int *crea_array_ent(int dil);
void guarda_resultados(int *arra, int nmpix,  FILE *fiiii);
FILE *apertura_archivo(int con);
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
int main(int argc, char **argv)
{
int i, j,  cont,  hit, no_hit, dimx, dimy, num_pix, factordim, NumAdq, N, NumPix, Niter, pixx, pixy;
int array_inp[3];
double  d, deltaX, deltaY, deltaP, betarnd, num_rand, mu, anchura_malla, delta, deltadis, deltaang, beta, deltaPaso;
clock_t comienzo, fin, co, cofi;
FILE  *fu1, *fprueba, *fu2;

#ifdef INPUT_TXT
fu1=fopen("input.txt", "r");

for(cont=0; cont<3; cont++)
{
fscanf(fu1, "%d", &array_inp[cont]);
printf("array_input[%d]=%d\n",cont, array_inp[cont]);


}
#endif
fclose(fu1);
factordim=2; 

N=array_inp[1];
NumAdq=array_inp[0]; /*Numero de angulos de adquisición*/
deltaang=2*PI/NumAdq;
anchura_malla=1; /*En unidades de medida, anchura del objeto*/
deltadis=0.001;
d=0.5;
deltaX=deltaY=0.01; /*Intervalo o paso del mesh grid*/
/*deltaP=0.01; Tamaño de cada detector*/
deltaPaso=0.001; /*Intervalo de propagación del rayo en cada iteración*/
delta=d; /*Distancia de la fuente al objeto*/
hit=0;
no_hit=0;
/*NumPix=(int)anchura_malla/deltaP; Numero de detectores*/
NumPix=array_inp[2];
deltaP=anchura_malla/NumPix;
dimx=(int)(anchura_malla/deltaX); /*Num de nodos en direccion x*/
dimy=(int)(anchura_malla/deltaY);
Niter=(int)(dimx*deltaX/deltaPaso);/*Numero de iteraciones/propagaciones de cada rayo*/

double *pos0 = malloc(2*sizeof(double));
double *pos1 = malloc(2*sizeof(double));
int *pix =malloc(2*sizeof(int));
int *resultados = malloc(NumPix*sizeof(unsigned long int)); /*Esta malloc...*/
 
srand(time(NULL));

/*Desde aquçi*/
//~ double *objeto = malloc(factordim*factordim*dimx*dimy*sizeof(double));
//~ objeto=crea_matriz(factordim*dimx, factordim*dimy, factordim*dimx/5, factordim*dimy/5);

//~ fu1=fopen("Matriz_objeto2.txt", "w+");
//~ guarda_array(objeto, factordim*dimx, factordim*dimy, fu1); /*Una vez guardado aqui se podría leer desde aquí*/
//~ free(objeto);
//~ fclose(fu1);
/*Hastaaqui*/

trataobj(fu1, factordim, 5, dimx, dimy, factordim*dimx/5, factordim/5*dimy); /*Crea objeto y lo guarda en un archivo de texto*/
fu1=fopen("Matriz_objeto2.txt", "r+");

double angmax=anchura_malla/(2*(anchura_malla+d)); 


for(cont=0; cont < NumAdq; ++cont) /*Iteracion de angulos */
{
co=clock();
beta=cont*deltaang;
resultados=crea_array_ent(NumPix); 
for(i=0; i<N; i++) /*Iteracion de todos los rayos para un mismo angulo*/
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
	pixx=-1;
	pixy=-1;



	for(j=0; j<Niter; j++) /*Recorrido de cada rayo*/

	{
	comienzo=clock();
	*(pos1)=*(pos0+1)*deltaPaso+*pos0;
	*(pos0)=*(pos1);
	pix=my_pixel(dimx, dimy, j, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0, beta); /*aquí se calcula un numero de pix coherente para distintos ángulos de adquisicion*/



    /*fseek(fu1, (pixx*dimy*factordim+pixy)*sizeof(double), SEEK_SET);  he quitado esto (pixx*dimy*factordim+pixy)*
	fscanf(fu1, "%lf", &mu);	Atenuacion lineal
	rewind(fu1);*/
	
if(pixx != *(pix) || pixy != *(pix+1)) /*Esto es un OR */
	{
	pixx=*(pix);
	pixy=*(pix+1);
	calcula_mu((pixx*dimy*factordim+pixy), &mu, fu1);
	num_rand=rand()/((double)RAND_MAX +1);

       /*printf("el angulo es %d la N es = %dla j es = %d y el pixel %d %d, cuya mu es %lf\n", cont, i,  j, *pix, *(pix+1), mu);*/




		if(num_rand > exp(-mu*deltaX)) /*Pasa o no el rayo*/
		{	no_hit=no_hit+1;
			printf("Rayo no llego\n");
			break;
	    }

	}

if(j== (Niter-1))
{
	num_pix=rint((*pos0-factordim*anchura_malla/2+anchura_malla/2)/deltaP); /*Num de detector al que llega el rayo*/
/*printf("numpix = %d, pos0 = %lf, pos0prim = %lf beta = %lf\n", num_pix, *pos0, *pos0-factordim*anchura_malla/2+anchura_malla/2, beta);*/

if(num_pix >= 0 && num_pix < NumPix)
/*printf("numpix es %d\n", num_pix);*/
	*(resultados+num_pix)=*(resultados+num_pix)+1;
printf("Rayo llego al pixel %d y el pixel tiene el valor de %d\n", num_pix, *(resultados+num_pix));
	hit=hit+1;
	}
	fin=clock();
	/*printf("Tiempo de ejecución de un rayo = %lf", ((double)(-comienzo+fin)/CLOCKS_PER_SEC));*/

}
}
printf("Final del ángulo %lf\n", beta);

fprueba=apertura_archivo(cont); /*Crea un archivo para cada ángulo*/

	guarda_resultados(resultados, NumPix, fprueba); /*Guarda los rayos detectados por cada pixel*/
	fclose(fprueba);
	free(resultados);
	 /*vuelvo a inicializar el array resultados*/
	 cofi=clock();
	 
	 //~ fu2=fopen("Cuenta de tiempo.txt", "a+");
	 
	 //~ fprintf(fu2, "Al angulo %d le ha costado %Lf\n", cont, (long double)(cofi-co));
	 //~ fclose(fu2);
}
printf( "\nHan llegado = %d \nNo han llegado %d de %lu rayos\n", hit, no_hit, N*NumAdq);





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
int *p = malloc(dilm*sizeof (*p));

	for(a=0; a<dilm; a++)
	*(p+a)=0;

	return p;

}
void trataobj(FILE *fui, int fac, int fac2, int dix, int diy, int deltax, int deltay)
{
	
int alm =fac*fac*dix*diy;
int a,b;
double objeto[alm];

fui=fopen("Matriz_objeto2.txt", "w+");

	

	for(a=0; a<dix; a++)
	{
		for(b=0; b<diy; b++)
		{		
			*(objeto+a*diy+b)=0.01;
		if(a < dix/2 + deltax && a > dix/2 -deltax && b < diy/2 +deltay && b > diy/2 - deltay) /*El objeto aquí está mal definido */
		*(objeto+a*diy+b)=5;
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
