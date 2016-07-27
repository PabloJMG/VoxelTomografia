#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#define deltaX 0.0005
#define deltaY 0.0005
#define deltaL 0.0005
#define deltaP 0.0001
#define NumPix 1000 /*Calcular numpix y numang automaticamente ¿?¿?*/
/*#define RANDINPUT*/
#define anchura_malla 1
#define INPUT_ARCHIVO

#define PI 3.141592653



double *import_data(double *datain);
int *crea_matriz(int dx, int dy);
unsigned long int *crea_array (int ded);
void guarda_matriz(int *res, int dix, int diy, FILE *fghq);
FILE *apertura_archivo(int co);
int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet);

int main(int argc, char **argv)
{

int dimx, dimy, i, cont, n, d, j, q,  numpixx2, numpixy2,  pasx, factordim, NumAng;

 int *numpix=malloc(2*sizeof(int));
double  yp, yl, division,  f, betarand, betamin, betamax, ang, beta;



double pos0rnd[2], pos1rnd[2];

NumAng=argv[1];
FILE *finp, *fout;
ang=2*PI/NumAng;
dimx=(int)anchura_malla/deltaX;
dimy=(int)anchura_malla/deltaY;
double deltaPaso;
deltaPaso=0.0001; /*Intervalo de propagacion de cada rayo*/
factordim=2;
pasx=(int)dimx*deltaX/deltaPaso;

int *matriz_resul = malloc(dimy*dimx*sizeof(*matriz_resul));
unsigned long int *datainp = malloc(NumPix*sizeof(*datainp));

matriz_resul=crea_matriz(factordim*dimx, factordim*dimy); /* Crea el mallado que va a guardar los resultados*/
datainp=crea_array(NumPix);

printf("Comienzo retroproyeccion\n");
getchar();
#ifdef RANDINPUT /*Inicializacion de un array con num enteros aleatorios*/
srand(time(NULL));
for(i=0; i<NumPix; i++)
{
datainp[i]=rand() % 11;

}
#endif

n=0;

#ifdef INPUT_ARCHIVO /*Coge los outputfiles del otro programa*/

for(cont=0; cont<NumAng; cont++)
{
beta=cont*ang;
finp=apertura_archivo(cont);
if(finp==NULL)
return 1;
printf("nuevo angulo\n");

for(i=0; i<NumPix; i++) /*Carga los datos de cada angulo de adq en datainp*/
{

fscanf(finp, "%lu", &datainp[i]);
printf("%lu\n", datainp[i]);

}




fclose(finp);


#endif


printf("Comienza el cuerpo de la simulación\n");
for(i=0; i<NumPix; i++)/*Cuerpo de toda la simulación, proyeccion de cada pixel */
{

if(datainp[i] != 0)
{
/*sabiendo la posicion del detector y la distancia lenslet-detector se tiene un cono angular dentro del cual se escoge un angulo aleatoriamente*/
yp=i*deltaP-deltaP/2;
division=yp/deltaL;
division=(int)division;
yl=division*deltaL+deltaL/2;
d=f=0.001;
betamin=(yp-yl-0.5*deltaP)/f;
betamax=(yp-yl+0.5*deltaP)/f;

for(j=0; j< datainp[i]; j++)
{
betarand=rand()/((double) RAND_MAX +1)*(betamax-betamin)+betamin;
pos0rnd[0]=yl;
pos0rnd[1]=betarand;
numpixx2=-1;
numpixy2=-1;

for(q=0; q<pasx; q++)
{
pos1rnd[0]=pos0rnd[0]+pos0rnd[1]*deltaPaso;
pos1rnd[1]=pos0rnd[1];
pos0rnd[0]=pos1rnd[0];
if(pos1rnd[0]<anchura_malla || pos1rnd[0]>0 )

	pos0rnd[0]=pos1rnd[0];


numpix=my_pixel(dimx, dimy, q, deltaPaso, deltaX, deltaY, factordim, anchura_malla, pos0rnd, beta); /*La funcion hace lo mismo en los dos programas, rota los pixeles en funcion del angulo de adquisicion*/
if(numpixx2!=*numpix && numpixy2!=*(numpix+1)) /*Este if probablemente sobre pero me libra de posibles abortos del programa por error */
{
/*Matriz_res= ¿?¿?¿? */
printf("El pixel es %d, %d,y el ang %lf\n", *numpix, *(numpix+1), beta);
*(matriz_resul+(*numpix*dimx)+*(numpix+1))=*(matriz_resul+(*numpix*dimx)+*(numpix+1))+1;
numpixx2=*numpix;
numpixy2=*(numpix+1);
}
}
}

}

}
/*Guarda todos los resultados en una matriz*/
fout=fopen("Resultados_Backprojection_Multiple files.txt", "w+");
guarda_matriz(matriz_resul, dimx, dimy, fout);
fclose(fout);
}
return 0;
}





int *crea_matriz(int dix, int diy)
{
int a, b;
int *r=malloc(dix*diy*sizeof(*r));
	for(a=0; a<dix; a++)

	{
		for(b=0; b<diy; b++)

		{

		*(r+a*dix+b)=0;


	}

		}
		return r;

}



double *import_data(double *datai) /*Al final no la uso o qué pasa ¿? */
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




}
FILE *apertura_archivo(int contador)
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "angulo_%d.txt", contador);

return fopen(nombrearchivo, "r");



}


int *my_pixel(int dix, int diy, int j, double deltaPas, double deltX, double deltY, int factodim, double anchura_mall, double *pos, double bet)
{
int *pix=malloc(2*sizeof(int));
int *pix2=malloc(2*sizeof(int));

		*pix=rint((j*deltaPas-anchura_mall/2)/deltX);
		*(pix+1)=rint(((*pos-anchura_mall/2)/deltY));
		*pix2=*pix*cos(bet)+*(pix+1)*sin(bet);
		*(pix2+1)=*(pix+1)*cos(bet)-*(pix)*sin(bet);
*pix=*pix2+factodim/2*dix;
*(pix+1)=*(pix2+1)+factodim/2*diy;

return pix;
}


