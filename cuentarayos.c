#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
int contador;
int suma;
int numero;
FILE *fin;
/*argv[2]=NumAng*/
int NumAng=20;
suma=0;
char nombrearchivo[25]="Resultados_Pixeles.txt";
for(contador=0; contador<NumAng; contador++)
{
sprintf(nombrearchivo, "angulo_%d.txt", contador);

fin=fopen(nombrearchivo, "r");

while(!feof(fin))
{

fscanf(fin, "%d", &numero);
suma=suma+numero;

}

}
printf("Han llegado %d rayos", suma);

}
