#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#define deltaX 0.02
#define deltaY 0.02
//~ #define NumPix 100 /*Calcular numpix y numang automaticamente ¿?¿?*/
/*#define RANDINPUT*/
//~ #define anchura_malla 1
#define INPUT_TXT
#define INPUT_ARCHIVO
#define PI 3.141592653


void guarda_matriz(int *res, int dix, int diy, FILE *fghq);
double *import_data(double *datain);
int *crea_matriz(int dx, int dy);
unsigned long int *crea_array (int ded);
int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet);
FILE *apertura_archivo(int co);
FILE *apertura_archivo2(int co, int co2, int co3);
int main(int argc, char **argv)
{

int dimx, dimy,  i, cont,  j, q,  numpixx2, numpixy2,  pasx, factordim, NumAng, conta, Num_Pix; /*Variables n, d no usadas*/
 //~ int *numpix=malloc(2*sizeof(int));
int *numpix=malloc(2*sizeof(int));
int array_inp[3];
double  yp, yl, division,  f, betarand, betamin, betamax, ang, beta, deltaPaso, deltaL, deltaP, anchura_malla, tiempo;
double pos0rnd[2], pos1rnd[2];
clock_t comienzo, fin;
FILE *finp, *fout, *fpru;
 #ifdef INPUT_TXT
finp=fopen("input.txt", "r");

for(cont=0; cont<3; cont++)
{
fscanf(finp, "%d", &array_inp[cont]);
printf("array_input[%d]=%d\n",cont, array_inp[cont]);


}
fclose(finp);
Num_Pix=array_inp[2];
NumAng=array_inp[0];

#endif 
anchura_malla=1;

conta=0;


ang=2*PI/NumAng;
dimx=(int)anchura_malla/deltaX;
dimy=(int)anchura_malla/deltaY;
deltaP=anchura_malla/Num_Pix;
deltaL=deltaP*5;
deltaPaso=0.0001; /*Intervalo de propagacion de cada rayo*/
factordim=2;
pasx=(int)dimx*deltaX/deltaPaso;


int *matriz_resul = malloc(dimy*dimx*sizeof(*matriz_resul));/*Esta función reserva 16 Mb de memoria*/
unsigned long int *datainp = malloc(Num_Pix*sizeof(*datainp)); 

matriz_resul=crea_matriz(factordim*dimx, factordim*dimy); /* Crea el mallado que va a guardar los resultados*/
datainp=crea_array(Num_Pix);

//~ printf("Comienzo retroproyeccion\n");
//~ getchar();
#ifdef RANDINPUT /*Inicializacion de un array con num enteros aleatorios*/
srand(time(NULL));
for(i=0; i<Num_Pix; i++)
{
datainp[i]=rand() % 11;

}
#endif

//~ n=0;

#ifdef INPUT_ARCHIVO /*Coge los outputfiles del otro programa*/

for(cont=0; cont<NumAng; cont++)
{
beta=cont*ang;
comienzo=clock();
finp=apertura_archivo(cont);
if(finp==NULL)
return 1;
printf("nuevo angulo\n");

for(i=0; i<Num_Pix; i++) /*Carga los datos de cada angulo de adq en datainp*/
{

fscanf(finp, "%lu", &datainp[i]);
printf("%lu\n", datainp[i]);

}



fclose(finp);

#endif


printf("Comienza el cuerpo de la simulación\n");
for(i=0; i<Num_Pix; i++)/*Cuerpo de toda la simulación, proyeccion de cada pixel */
{

if(datainp[i] != 0)
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


for(j=0; j< datainp[i]; j++)
{
 betarand=rand()/((double) RAND_MAX +1)*(betamax-betamin)+betamin; /*Esto no debería ser de signo contrario ¿?*/
/*betarand=0.05;*/
pos0rnd[0]=yl;
pos0rnd[1]=betarand;
numpixx2=-1;
numpixy2=-1;
//~ fpru=apertura_archivo2(i, j, cont);
//~ fprintf(fpru,  "%lf   %lf	%d		%d\n", pos0rnd[0], pos0rnd[1], numpixx2, numpixy2);
for(q=0; q<pasx; q++)
{
pos1rnd[0]=pos0rnd[0]+pos0rnd[1]*deltaPaso;
pos1rnd[1]=pos0rnd[1];
pos0rnd[0]=pos1rnd[0];
//~ if(pos1rnd[0]<anchura_malla || pos1rnd[0]>0 ) /*Este if para qué ¿?*/

	//~ pos0rnd[0]=pos1rnd[0];


numpix=my_pixel(dimx, dimy, q, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0rnd, beta); 


/*printf("Pos0 = %lf, angulo = %lf, rango de angulos = %lf, %lf, pix= %d, piy=%d\n", pos0rnd[0], pos0rnd[1], betamax, betamin, numpix[0], numpix[1]);*/
if(numpixx2!=*numpix || numpixy2!=*(numpix+1)) /*Este if probablemente sobre pero me libra de posibles abortos del programa por error */
{
/*Matriz_res= ¿?¿?¿? */
/*printf("El pixel es %d, %d,y el ang %lf\n", *numpix, *(numpix+1), beta);*/

//~ *(matriz_resul+(*numpix)*dimx+*(numpix+1))+=1;
numpixx2=*numpix;
numpixy2=*(numpix+1);
 //~ *(matriz_resul+numpixx2*dimx+numpixy2)+=1;
 *(matriz_resul+numpixy2*dimx*factordim+numpixx2)+=1;
conta+=1;
//~ fprintf(fpru, "%lf   %lf	%d		%d\n", pos0rnd[0], pos1rnd[1], numpixx2, numpixy2);
}
}
}
//~ fclose(fpru);
}

}
fin=clock();
tiempo=difftime(fin, comienzo);
fpru=fopen("Calculo_tiempo.txt", "a+");

fprintf(fpru, "%lf\n", tiempo);
fclose(fpru);

/*Guarda todos los resultados en una matriz*/

}
char archivo[50]="Resultados_BP.txt";
fout=fopen(archivo, "w+");


guarda_matriz(matriz_resul,factordim*dimx, factordim*dimy, fout);
fclose(fout);
FILE *pipeplot=popen("gnuplot -persist", "w");
/*fprintf(pipeplot, "set title /lo que sea/"   */
fprintf(pipeplot, "set pm3d map\n");
fprintf(pipeplot, "splot '/home/pablo/Resultados_C/FP_BP/Resultados_BP.txt' matrix\n");
fflush(pipeplot);
sleep(5);
fprintf(pipeplot, "quit\n");
pclose(pipeplot);
//~ printf("conta=%d\n", conta);
//~ getchar();
return 0;
}





int *crea_matriz(int dix, int diy)
{
int a;

int dim=dix*diy;
 int *r=malloc(dim*sizeof(int));

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
		free(r);

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
	  unsigned long int *p = malloc(di*sizeof(*p));

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
printf("Final de la subrutina guarda_matriz");



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

int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet)
{
static int pix[2];
int pix2[2];

		*pix=rint((j*deltaPas-anchura_mall/2)/deltX); 
		*(pix+1)=rint(((*pos-anchura_mall/2)/deltY));
		*pix2=*pix*cos(bet)-*(pix+1)*sin(bet);
		*(pix2+1)=*(pix+1)*cos(bet)+*(pix)*sin(bet);
*pix=*pix2+factodim/2*dix;
*(pix+1)=*(pix2+1)+factodim/2*diy;

return pix;
}


