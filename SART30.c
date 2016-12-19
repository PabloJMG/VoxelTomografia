
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.141592
#define DOSPI
#define COEF_COM

void *carga_archivo(int NumAng, int NumPix, float *b);

FILE *apertura_archivo(int contador);

int main(int argc, char **argv)
{

int NumPix, cont, cont2, cont3, cont4,cont5, cont6,dimx, dimy, pixx, pixy, NumAng, i, valu, dimtotal, NumRayos;
float beta, cosbeta, factora, sinbeta, deltaang, lambda, deltax, deltay, deltaf, posx, posy, posx2, posy2, val, diferencia, cant;
double angulos[10];
FILE *finp;
clock_t inicio, final;


factora=20;
dimx=dimy=32;
int aint[3];
lambda=0.1; /*Parámetro de corrección*/
finp=fopen("inputint.txt", "r");
i=0;
for(cont=0; cont < 3; cont++)
	{		
	fscanf(finp, "%d", &aint[cont]);
	printf("%d\n", aint[cont]);
	}
fclose(finp);

#ifdef DOSPI
NumAng=aint[0];
#endif


#ifdef ABANICO /*Lo voy a dejar en version naif*/

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

NumPix=aint[2];

dimtotal=dimx*dimy;
NumRayos=NumAng*NumPix;


float *a=calloc(dimx*dimy*NumPix*NumAng, sizeof(float));
float *b = calloc(NumPix*NumAng, sizeof(float)); /*cambiado hoy */
 float *x=calloc(dimx*dimy, sizeof(float));
 float *xantes=calloc(dimx*dimy, sizeof(float));
 float *betar=calloc(NumPix*NumAng, sizeof(float)); /*Sumatorio de los coeficientes a_ij a lo largo de toda la imagen referidos al mismo pixel del mismo ángulo */
 float *gamma=calloc(dimtotal*NumAng, sizeof(float)); /* Sumatorio de los coeficientes a_ij   */
 
  deltaang=2*pi/NumAng;
  //~ lambda=2;
  /* Inicialización de la imagen x */
  
  for(cont=0; cont<NumAng; cont++)
  {
	  	 *(x+cont)=0;
  }
  
  for(cont=0; cont<(dimtotal*NumAng); cont++)
  {
	  
	  *(gamma+cont)=1;
  }

#ifdef COEF_INT
for(cont=0; cont<NumAng; cont++)
{	
	 #ifdef DOSPI 
	 beta=cont*deltaang;
	 #endif
	 
	  #ifdef ABANICO
	 beta=angulos[cont];
	 #endif
		 cosbeta=cos(beta);
	 sinbeta=sin(beta);
	for(cont2=0; cont2<NumPix; cont2++)
	{
		//~ finp=apertura_archivo2(cont, cont2);
		for(cont3=0; cont3< dimx; cont3++)
		{
			  posx=cont3;	
			 posy=cont2*dimy/NumPix; /*Error aqui*/
			 //~ printf("posyprevia=%f, cont2=%d, dimy=%d, NumPix=%d\n",posy, cont2, dimy,NumPix);
			 posx2=posx-dimx/2;
			 posy2=posy-dimy/2;
			 posx=posx2*cosbeta + posy2*sinbeta; /*inventada a toda prisa*/
			 posy=posy2*cosbeta - posx2*sinbeta;
			 posx=posx+dimx/2;
			 posy=posy+dimy/2;
			 pixx=trunc(posx);
			 pixy=trunc(posy);
			 *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+pixy+cont3*dimy)+=1;
			
			 for(cont4=0; cont4<dimy; cont4++)
			 {
			 //~ fprintf(finp, "%f ", *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+cont4+cont3*dimy));
			 
		 }
		 //~ fprintf(finp,"\n");
		}
		
		//~ fclose(finp);
	}
	
	
}

#endif


#ifdef COEF_COM
 for(cont=0; cont<NumAng; cont++) /*inicializar la matriz de coeficientes a*/
 { 		
	 #ifdef DOSPI 
	 beta=cont*deltaang;
	 #endif
	 
	 #ifdef ABANICO
	 beta=angulos[cont];
	 #endif
	 cosbeta=cos(beta);
	 sinbeta=sin(beta);
	 for(cont2=0; cont2<NumPix; cont2++)
	{ //~ {finp=apertura_archivo2(cont, cont2);
		 for(cont3=0; cont3<(dimx-1); cont3++) /*El limite del contador puesto así para que los píxeles limítrofes tengan menos peso*/
		 {	
			 //~ for(cont4=0; cont4<(dimy-1); cont4++)
			 //~ {
			  posx=cont3;	
			 posy=cont2*dimy/NumPix; /*Error aqui*/
			 //~ printf("posyprevia=%f, cont2=%d, dimy=%d, NumPix=%d\n",posy, cont2, dimy,NumPix);
			 posx2=posx-dimx/2;
			 posy2=posy-dimy/2;
			 posx=posx2*cosbeta + posy2*sinbeta; /*inventada a toda prisa*/
			 posy=posy2*cosbeta - posx2*sinbeta;
			 posx=posx+dimx/2;
			 posy=posy+dimy/2;
			 pixx=trunc(posx);
			 pixy=trunc(posy);
			 deltax=posx-pixx;
			 deltay=posy-pixy;
	
			 *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+pixy+cont3*dimy)+=(2-deltax-deltay)/factora;/*Pensar esta mierda---- mal indizado*/
			 *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+(pixy+1)+cont3*dimy)+=(1-deltax+deltay)/factora;
			 *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+pixy+(cont3+1)*dimy)+=(1-deltay+deltax)/factora;
			 *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+(pixy+1)+(cont3+1)*dimy)+=(deltax+deltay)/factora;/*Dará SEGFAULT segurp*/
			 		//~ printf("posx=%f, posy=%f, pixx = %d, pixy=%d", posx, posy, pixx, pixy);
			
			//~ printf("a=%lf", *(a+cont*NumPix*dimtotal+cont2*dimtotal+cont3*dimy);
			 //~ for(cont4=0; cont4<dimy; cont4++)
			 //~ {
			//~ printf( "%lf\n", *(a+cont*NumPix*dimx*dimy+cont2*dimx*dimy+cont4+cont3*dimy));
			 
		 //~ }
		 //~ fprintf(finp,"\n");
		 //~ getchar();
	 }
	
	//~ fclose(finp);
	
	 }
	 
	 
 }
 #endif
/*Carga los resultados de las proyecciones desde los archivos a un array*/
 carga_archivo(NumAng, NumPix, b);
/* Se calculan las cantidades beta y gamma */
inicio=clock();
for(cont=0; cont<NumAng; cont++)
 {
for(cont2=0; cont2<NumPix; cont2++)
{
	for(cont3=0; cont3<dimx; cont3++)
	{
		for(cont4=0; cont4<dimy; cont4++)
		{
			
			*(betar+cont*NumPix+cont2)+=*(a+cont*NumPix*dimtotal+cont2*dimtotal+cont3*dimy+cont4);
		
			
		}
	}

//~ printf("beta=%lf\n", *(betar+cont*NumPix+cont2));
}	
	 
	 
 }
 
 for(cont=0; cont<NumAng; cont++)
	{
		for(cont2=0; cont2<dimx; cont2++)
		{
			for(cont3=0; cont3<dimy; cont3++)
			
			{
				for(cont4=0; cont4<NumPix; cont4++)
				{
					*(gamma+cont*dimtotal+cont2*dimy+cont3)+=*(a+cont*NumPix*dimtotal+cont4*dimtotal+cont2*dimy+cont3); /*Esta bien esta inicializacion¿? */
					
					
				}
			//~ printf("gamma=%lf, cont=%d, cont2=%d, cont3=%d\n", *(gamma+cont*dimtotal+cont2*dimy+cont3), cont, cont2, cont3 );
			}
			
		}
		
	}
//~ getchar();
i=0;
deltaf=1;
  while(fabs(deltaf)>0.0001)
 {
 deltaf=0;
	 for(cont=0; cont<NumAng; cont++)
	 {
		for(cont2=0; cont2<dimx; cont2++)
		{
			for(cont3=0; cont3<dimy; cont3++)
			{
				for(cont4=0; cont4 < NumPix; cont4++)
				{
					cant=0;
					
					for(cont5=0; cont5<dimx; cont5++)
					{
						
						for(cont6=0; cont6<dimy; cont6++)
						
								{
								cant+=*(a+cont*NumPix*dimtotal+cont4*dimtotal+cont5*dimy+cont6)*(*(x+cont5*dimy+cont6));
								//~ printf("cant=%lf\n", cant);
								}		
					}
				
				diferencia+=(*(b+cont*NumPix+cont4)-cant)*(*(a+cont*NumPix*dimtotal+cont4*dimtotal+cont2*dimy+cont3))/(*(betar+cont*NumPix+cont4));
			//~ printf("diferencia=%lf\n", diferencia);
			}
			//~ getchar();
			
		*(x+cont2*dimy+cont3)+=lambda*diferencia/(*(gamma+cont*dimtotal+cont2*dimy+cont3));
		
		deltaf+=diferencia/(*(gamma+cont*dimtotal+cont2*dimy+cont3));
		//~ printf("deltaf=%lf, difer=%lf, gamma=%lf cont =%d, cont2=%d\n", deltaf, diferencia, *(gamma+cont*dimtotal+cont2*dimy+cont3), cont, cont2);
		//~ getchar();
		diferencia=0;
		}	
		
	}

	 
 }
 
 //~ deltaf=deltaf*deltaf; /* arg_{min} (b-Ax)²  */
 	i+=1;
	 printf("Deltaf=%lf y n de iteraciones =%d\n", deltaf, i);
	 //~ getchar();
 
}
final=clock();
printf("Tardo %lf segundos para %d angulos, %d pixeles y tamaño %d x %d\n", (double)(final-inicio)/CLOCKS_PER_SEC, NumAng, NumPix, dimx, dimy);
finp=fopen("Resultados_SART.txt", "w+");

for(cont=0; cont<dimy; cont++)
{
	
	for(cont2=0; cont2<dimx; cont2++)
	{
		
		
		fprintf(finp, "%f   ", *(x+dimx*cont+cont2)); 
	}
	fprintf(finp, "\n");
}
//~ *(g+m)+=*(a+m*dimI+i)*(*(p+i)-(*(at+m*dimI+i)*(*(g+ m*dimI+i)))/((*(at+m*dimI+i))*(*(a+m*dimI+i)));
fclose(finp);


return 0;
}


void *carga_archivo(int NumAng, int NumPix, float *b)
{
	int cont, cont2, valu;
	
	FILE *finp;
	
 for(cont=0; cont<NumAng; cont++) /*  Cargar b (esto está mal)*/
 {
	 
	 finp=apertura_archivo(cont);
	 
	 for(cont2=0; cont2<NumPix; cont2++)
	 {
	fscanf(finp, "%d", &valu);
	printf("val=%d\n", valu);
	*(b+cont*NumPix+cont2)=log(1+valu); /*Esto aquí lo he cambiado*/
	//~ *(b+cont*NumPix+cont2)=valu;
	printf("val=%d, log(val)=%f\n", valu, *(b+cont*NumPix+cont2));
	
	 }
	 
	 fclose(finp);
 }
	
	//~ printf("Archivo cargado, presione alguna tecla\n");
	//~ getchar();
}

FILE *apertura_archivo(int contador)
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "angulo_%d.txt", contador);

if(fopen(nombrearchivo, "r")==NULL)
{
	perror("El archivo no existe o no se puede abrir\n");
}

return fopen(nombrearchivo, "r");



}
