


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
int *crea_array_ent(int dil);
void guarda_resultados(int *arra, int nmpix, FILE *fiiii);
/*double *matrixmul(double aaa, double bbb);*/
/*#define anchura_malla 1
#define deltaY 0.001
#define deltaX 0.001
#define dimx anchura_malla/deltaX
#define dimy anchura_malla/deltaY*/
int main(int argc, char **argv)
{
int i, j, N,NumPix, hit, no_hit, dimx, dimy, pixx, pixy, num_pix;
double  d, deltaX, deltaY, deltaP, betarnd, num_rand, mu, anchura_malla, delta;

anchura_malla=1;
FILE *fout, *fpix, *fu1, *fprueba;
d=0.5;
N=1000000; /*Número de iteraciones*/
deltaX=deltaY=0.01;
deltaP=0.001;
delta=d;
hit=0;
no_hit=0;
NumPix=(int)anchura_malla/deltaP;
dimx=(int)(anchura_malla/deltaX);
dimy=(int)(anchura_malla/deltaY);

double *pos0 = malloc(2*sizeof(double));

double *pos1 = malloc(2*sizeof(double));

double *objeto = malloc(dimx*dimy*sizeof(double));
objeto=crea_matriz(dimx, dimy, dimx/10, dimy/10);
/*objeto=crea_matriz(dimx, dimy, dimx/10, dimy/10); Este objeto se define muy mal, mirar la función que lo genera */
int *resultados = malloc(NumPix*sizeof(int));
resultados=crea_array_ent(NumPix); /*Otra funcion crea_array con retorno de entero*/
double angmax=anchura_malla/(2*(anchura_malla+d));
fout=fopen("resultados_codeblocks.txt", "w+");
fpix=fopen("pixeles.txt", "w+");
srand(time(NULL));
fprueba=fopen("Resultados_Pixeles2.txt", "w+");
printf("Comienzo proyeccion\n");
getchar();

for(i=0; i<dimx; i++)
{for(j=0; j<dimy; j++)
	{

	printf("%lf\n",*(objeto +i*dimx+j));
	/*printf("Da el fallo en el elemento %d, %d; el elemento es: %lf \n", j, i,*(objeto +i*dimx+j));*/
	}

}

fu1=fopen("Matriz_objeto.txt", "w+");

/*Quitado el pos0[0] > anchura_malla pq ya en la def de angmax impongo esta cond */
for(i=0; i<N; ++i) /*Proyecciones de todos los rayos*/
{
	num_rand=rand()/((double)RAND_MAX +1); /*Desbordamiento entero aquí*/
	betarnd=-angmax+2*angmax*num_rand;

	*pos0 = anchura_malla/2;
	 *(pos0+1)= betarnd;

	 guarda_array(objeto, dimx, dimy, fu1);

	*pos1=*(pos0+1)*delta+*pos0;
	*pos0=*pos1;
	*(pos1+1)=*(pos0 +1 );
fprintf(fout, "%f %f %f %d %d\n", pos0[0], pos0[1],betarnd, i, 0);


	for(j=0; j<dimx;++j)

	{

	*pos1=*(pos0+1)*deltaX+*pos0;
	*pos0=*pos1;

		pixx=dimx-j-1;
		pixy=dimy-rint((*pos1/deltaX));


		mu=*(objeto + dimx*pixx + pixy);
		num_rand=rand()/((double)RAND_MAX +1);
		printf("mu=%lf, numrand=%lf\n", mu, num_rand);
        printf("la j es = %d y el pixel %d %d", j, pixx, pixy);


	fprintf(fpix, "%d %d\n", pixx, pixy);

		if(num_rand > exp(-mu*deltaX))
		{no_hit=no_hit+1;

			break;
	}

	fprintf(fout, "%f %f %f %d %d\n", *pos0, *(pos0+1),betarnd, i, j);
}
if(pixx== 0)
{
	num_pix=rint(*pos0/deltaP);
	resultados[num_pix]++;

	hit=hit+1;
	}
}

printf( "\nHan llegado = %d \nNo han llegado %d de %d rayos\n", hit, no_hit, N);
	guarda_resultados(resultados, NumPix, fprueba);
	fclose(fprueba);
	fclose(fpix);
	fclose(fout);
fclose(fu1);

	printf("Va bene\n");
	getchar();

	return 0;
}


double *crea_matriz(int dix, int diy, int deltax, int deltay)
{
int a, b;
double *r=malloc(dix*diy*sizeof(*r));
	for(a=0; a<dix; a++)

	{
		for(b=0; b<diy; b++)

		{
		if(a < dix/2 + deltax && a > dix/2 -deltax && b < diy/2 +deltay && b > diy/2 - deltay) /*El objeto aquí está mal definido */
		*(r+a*dix+b)=10000;
		else
		*(r+a*dix+b)=0;

	}

		}
		return r;

}
double *crea_array(int di)
{

	int a;
    double *p =malloc(di*sizeof *p);

	for(a=0; a<di; a++)
	{
	p[a]=0;
}
	return p;

}

int *crea_array_ent(int dilm)
{

	int a;
	  int *p = malloc(dilm*sizeof *p);

	for(a=0; a<dilm; a++)
	*(p+a)=0;

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
 fprintf(fau, "%lf", r);


}

fprintf(fau, "\n");
}

/*fclose(fau);*/
}

void guarda_resultados(int *arrayn, int npi, FILE* finp)
{
int cont;
int res;

for(cont=0; cont<npi; cont++)
{
res=*(arrayn + cont);
fprintf(finp, "%d\n",  res);




}


}
