/*
*
*
*
*
*
*
*
*
*
*
*
*
*
*/

/*En MatLab definía un vector con números aleatorios enteros que introducía como input */


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
#define NumPix 1000
/*#define RANDINPUT*/
#define anchura_malla 1
#define INPUT_ARCHIVO
double *import_data(double *datain);
int *crea_matriz(int dx, int dy);
int *crea_array (int ded);
void guarda_matriz(int *res, int dix, int diy, FILE *fghq);
int main(int argc, char **argv)
{

int dimx, dimy, i, j, q, n, numpixx2, numpixy2, numpixy,numpixx,ch;



double  yp, yl, division, d, f, betarand, betamin, betamax;

char file_name[250]="Resultados_Pixeles2.txt";

double pos0rnd[2], pos1rnd[2];
size_t result;
int *buffer;
FILE *finp, *fout;

dimx=(int)anchura_malla/deltaX;
dimy=(int)anchura_malla/deltaY;
int *matriz_resul = malloc(dimy*dimx*sizeof(*matriz_resul));
int *datainp = malloc(NumPix*sizeof(*datainp));

matriz_resul=crea_matriz(dimx, dimy); /* Crea el mallado que va a guardar los resultados*/
datainp=crea_array(NumPix);

printf("Comienzo retroproyeccion\n");
getchar();
#ifdef RANDINPUT
srand(time(NULL));
for(i=0; i<NumPix; i++)
{
datainp[i]=rand() % 11;

}
#endif

n=0;

#ifdef INPUT_ARCHIVO
printf("Por favor, introduzca el nombre del archivo de input.\n");
gets(file_name);

finp=fopen(file_name, "r");


if(finp == NULL)
{
printf("El nombre del archivo es incorrecto o no existe\n");
exit(EXIT_FAILURE);

}

/*Calcular tamaño de archivo input
fseek(finp, 0, SEEK_END);
n=ftell(finp);
fseek(finp,0, SEEK_SET);

printf("Hay %d datos en el archivo.\n", n);
rewind(finp);


buffer=(int*)malloc(NumPix*sizeof(int*));
result=fread(buffer, sizeof(buffer[0]), n, finp); Definir  result ¿?¿?*/





for(i=0; i<NumPix; i++)
{
/*datainp[i]=getc(finp);*/ /*Fallo aquí seguro */
fscanf(finp, "%d", &datainp[i]);
printf("%d\n", datainp[i]);

}


fclose(finp);

/*Calcular la longitud del vector antes y luego que lo guarde donde sea ¿?¿? */
#endif

printf("Comienza el cuerpo de la simulación\n");
for(i=0; i<NumPix; i++)/*Cuerpo de toda la simulación */
{

if(datainp[i] != 0)
{
yp=i*deltaP-deltaP/2;
division=yp/deltaL;
division=(int)division; /*Probablemente la función aquí sea rint */
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

for(q=0; q<dimx; q++)
{
pos1rnd[0]=pos0rnd[0]+pos0rnd[1]*deltaX;
pos1rnd[1]=pos0rnd[1];
pos0rnd[0]=pos1rnd[0];
if(pos1rnd[0]<anchura_malla || pos1rnd[0]>0 )

	pos0rnd[0]=pos1rnd[0];
numpixx=q;
numpixy=(int)(pos0rnd[0]/deltaY);

if(numpixx2!=numpixx && numpixy2!=numpixy) /*Este if tiene la condición bien definida */
{
/*Matriz_res= ¿?¿?¿? */
*(matriz_resul+numpixx*dimx+numpixy)=*(matriz_resul+numpixx*dimx+numpixy)+1;
numpixx2=numpixx;
numpixy2=numpixy;
}
}
}

}

}
fout=fopen("Resultados_Backprojection.txt", "w+");
guarda_matriz(matriz_resul, dimx, dimy, fout);
fclose(fout);

return 0;
}



/*Contar las llaves */


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

int *crea_array(int di)
{

	int a;
	  int *p = malloc(di*sizeof(*p));

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

