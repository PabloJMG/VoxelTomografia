
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/types.h>
#include <unistd.h>
#define pi 3.141592
#define DOSPI
#define COEF_COM

void *carga_archivo(int NumAng, int NumPix, float **b);

FILE *apertura_archivo(int contador);
FILE *apertura_archivo2(int cont1, int cont2);
void saltar_com(FILE *fin);
struct nonzer{
	
	int nonnum;
	int **pixeles;
	float *as; /*Necesario esto --- se podria liberar toda la matriz a*/
	
};
#ifdef COEF_COM
struct coef{
	
	float **a;
	
};
#endif

#ifdef COEF_INT
struct coef{
	
int  **a;
	
};
#endif
/*Vas a tener que cambiar todo por el factordim  ---- DONE*/ 

/* Cambiar el modo de guardar los pixeles distintos de cero, o al menos revisarlo*/
int main(int argc, char **argv)
{

int NumPix, cont, cont2, cont3, cont4, cont5, cont6, dimx, dimy, pixx, pixy, NumAng, i, valu, dimtotal, NumRayos, nth, limx, limy, con2, con3, pox, poy;
float beta, cosbeta, factora, sinbeta, deltaang, lambda, deltax, deltay, deltaf, posx, posy, posx2, posy2, val, diferencia, cant, sumpix, factorpos;
double angulos[10];
char filename[100];
FILE *finp;
clock_t inicio, final;

int factordim=2;
factora=0.02;
//~ dimx=dimy=64;
finp=fopen("input/dimSART.txt", "r");
fscanf(finp, "%d", &dimx);
fscanf(finp, "%d", &dimy);
fclose(finp);
int aint[3];
lambda=0.5; /*Parámetro de corrección*/
finp=fopen("input/inputint.txt", "r");
i=0;

for(cont=0; cont < 3; cont++)
	{
		saltar_com(finp);		
	fscanf(finp, "%d", &aint[cont]);
	printf("%d\n", aint[cont]);
saltar_com(finp);
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
	//~ printf("%lf\n", angulos[cont]);
	
}
fclose(finp);
//~ getchar();
#endif

NumPix=aint[2];
dimtotal=factordim*factordim*dimx*dimy;
NumRayos=NumAng*NumPix;
limx=factordim*dimx;
limy=factordim*dimy;
//~ struct nonzer{
	
	//~ int nonnum;
	//~ int pixeles[2][NumRayos];
	//~ float as[NumRayos]; /*Necesario esto --- se podria liberar toda la matriz a*/
	
//~ };
struct coef coe[NumAng][NumPix]; /*Necesario guardar todas las matrices de coeficientes o con el struct de abajo basta ¿? */
struct nonzer nonum[NumAng][NumPix];


for(cont=0; cont<NumAng; cont++)
{
	
	for(cont2=0; cont2<NumPix; cont2++)
	{
		
		coe[cont][cont2].a=malloc(limx*sizeof(float));
		for(cont3=0; cont3<limx; cont3++)
		{
			
			coe[cont][cont2].a[cont3]=malloc(limy*sizeof(float));
			for(cont4=0; cont4<limy; cont4++)
			{
				coe[cont][cont2].a[cont3][cont4]=0;
				
			}
		}
		
	}
	
}
factorpos=(float)dimy/NumPix;  /* factordim incluido, necesario ¿?*/
printf("factorpos=%f\n", factorpos);
//~ getchar();
		 //~ float *a=calloc(factordim*factordim*dimx*dimy*NumPix*NumAng, sizeof(float));
		 //~ float *b = callocc(NumPix*NumAng, sizeof(float)); /*cambiado hoy */
		 
float **b, **betar, **x;
struct coef gamma[NumAng];
b=calloc(NumAng, sizeof(float));
//~ gamma=calloc(NumAng, sizeof(float));
betar=calloc(NumAng, sizeof(float));
x=calloc(limx, sizeof(x));
				for(cont=0; cont<NumAng; cont++){
							b[cont]=calloc(NumPix, sizeof(*b));		
							betar[cont]=calloc(NumPix, sizeof(*betar));
							gamma[cont].a=calloc(limx, sizeof(*gamma));
									for(cont2=0; cont2<limx; cont2++) {
										
										gamma[cont].a[cont2]=calloc(limy, sizeof(*gamma[cont].a));	/* Sumatorio de los coeficientes a_ij   */
										x[cont2]=calloc(limy, sizeof(*x));
										
										}
				} 

  deltaang=2*pi/NumAng;


/* Alocacion de espacio para pixeles.nonum y as.nonum  */

		
		for(cont=0; cont<NumAng; cont++)
		{
			
			for(cont2=0; cont2<NumPix; cont2++)
			
			{
				 nonum[cont][cont2].as=malloc(dimtotal*sizeof(float));
				 nonum[cont][cont2].pixeles=malloc(2*sizeof(int *));
				 for(cont3=0; cont3<2; cont3++)
				 
				 nonum[cont][cont2].pixeles[cont3]=malloc(dimtotal*sizeof(int));
				
			}
		}
#ifdef COEF_INT /*No tocada mucho, no es la que suelo usar */
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
		finp=apertura_archivo2(cont, cont2);
		nth=0;
		for(cont3=0; cont3< dimx; cont3++)
		{
			  //~ posx=dimx-cont3;	
			  posx=factordim*dimx-(factordim-1)*dimx/2 -cont3;
			 posy=cont2*factorpos+(factordim-1)*dimy/2; /*Error aqui*/
			 //~ printf("posyprevia=%f, cont2=%d, dimy=%d, NumPix=%d\n",posy, cont2, dimy,NumPix);
			 posx2=posx-factordim*dimx/2;
			 posy2=posy-factordim*dimy/2;
			 posx=posx2*cosbeta + posy2*sinbeta; /*inventada a toda prisa*/
			 posy=posy2*cosbeta - posx2*sinbeta;
			 posx=posx+factordim*dimx/2;
			 posy=posy+factordim*dimy/2;
			 pixx=trunc(posx);
			 pixy=trunc(posy);
			
			coe[cont][cont2].a[pixx][pixy]+=1;
			 for(cont4=0; cont4<dimy; cont4++)
			 {
			
							 if(coe[cont][cont2].a[cont3][cont4]!=0) /* Se tiene que cambiar */
							 {

								 	nonum[cont][cont2].pixeles[0][nth]=pixx;
									nonum[cont][cont2].pixeles[1][nth]=pixy;
									nonum[cont][cont2].as[nth]=coe[cont][cont2].a[cont3][cont4];
									nth+=1;
									fprintf(finp, "%d %d %lf\n", cont3, cont4, coe[cont][cont2][cont3][cont4]);
							 }
		 }
		 //~ getchar();
		 		 //~ nonum[cont][cont2].nonnum= nth;
		 		 //~ printf("nth=%d\n", nth);
		 fprintf(finp,"\n");
		 //~ getchar();
		}
		 nonum[cont][cont2].nonnum= nth;
		 		 //~ printf("nth=%d\n", nth);
		 		 //~ getchar();
		fclose(finp);
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
	 for(cont2=0; cont2<NumPix; cont2++){
			//~ finp=apertura_archivo2(cont, cont2);
			nth=0;
						 for(cont3=0; cont3<dimx; cont3++) /*El limite del contador puesto así para que los píxeles limítrofes tengan menos peso*/
						 {	
							 //~ for(cont4=0; cont4<(dimy-1); cont4++)
							 //~ {
							  posx=factordim*dimx-(factordim-1)*dimx/2-cont3;	
							 posy=factordim*dimy-(factordim-1)*dimy/2 -cont2*factorpos; /*Error aqui*/

							 posx2=posx-factordim*dimx/2;
							 posy2=posy-factordim*dimy/2;
							 posx=posx2*cosbeta + posy2*sinbeta; /*inventada a toda prisa*/
							 posy=posy2*cosbeta - posx2*sinbeta;
							 posx+=factordim*dimx/2;
							 posy+=factordim*dimy/2;
							 //~ posy=2.0;
							 pixx=trunc(posx);
							 pixy=trunc(posy);

							 deltax=posx-pixx;
							 deltay=posy-pixy;
					
						
							 coe[cont][cont2].a[pixx][pixy]+=(2-deltax-deltay)/factora;
							 coe[cont][cont2].a[pixx][pixy+1]+=(1-deltax+deltay)/factora;
							 coe[cont][cont2].a[pixx+1][pixy]+=(1-deltay+deltax)/factora;
							 coe[cont][cont2].a[pixx+1][pixy+1]+=(deltax+deltay)/factora;
							 		//~ printf("posx=%f, posy=%f, pixx = %d, pixy=%d", posx, posy, pixx, pixy);
							
						
							 for(cont4=0; cont4<dimy; cont4++)
							 {
							//~ fprintf(finp,  "%d %d %lf\n", cont3, cont4, coe[cont][cont2].a[cont3][cont4]);
							pox=rint(factordim*dimx-(factordim-1)*dimx/2-cont3);
							poy=rint((factordim-1)*dimy/2+cont4);
							if(coe[cont][cont2].a[pox][poy]>0.00001) /* ¿?¿?¿?*/
							{
										nonum[cont][cont2].pixeles[0][nth]=pox;
										nonum[cont][cont2].pixeles[1][nth]=poy;
										nonum[cont][cont2].as[nth]=coe[cont][cont2].a[pox][poy];
										nth+=1;
							
							//~ printf("NumAng =%d, numpix = %d, x =%d, pixx = %d, y=%d,  pixy =%d, coe =%f \n", cont, cont2, cont3, pox, cont4, poy, coe[cont][cont2].a[pox][poy]);
							//~ getchar();	
							}
							 
						 }
						 //~ fprintf(finp,"\n");
						 //~ getchar();
						 
					 }
					 nonum[cont][cont2].nonnum= nth;		
					 printf("nth=%d\n", nth);
	//~ fclose(finp);
	
	 }
	 
	 
 }
 #endif
 getchar();
 for(cont=0; cont<NumAng; cont++)
 {
	 for(cont2=0; cont2<NumPix; cont2++)
	 { finp=apertura_archivo2(cont, cont2);
		 for(cont3=0; cont3<limx; cont3++)
		 {
			 for(cont4=0; cont4 < limy; cont4++)
			 {
				 
				 fprintf(finp, "%d %d %f\n", cont3, cont4, coe[cont][cont2].a[cont3][cont4]);
				 
			 }
		//~ fprintf(finp, "\n");
		
		 }
		 fclose(finp);
	 }
 }
/*Carga los resultados de las proyecciones desde los archivos a un array*/
 carga_archivo(NumAng, NumPix, b);
/* Se calculan las cantidades beta y gamma */
//~ printf("comienza iteración, presionar tecla\n");
//~ getchar();
inicio=clock();
for(cont=0; cont<NumAng; cont++)
 {
		for(cont2=0; cont2<NumPix; cont2++)
		{
			for(cont3=(factordim-1)*dimx/2; cont3<(factordim+1)*dimx/2; cont3++)
			{
				for(cont4=(factordim-1)*dimy/2; cont4<(factordim+1)*dimy/2; cont4++)
				{
				
							betar[cont][cont2]+=(coe[cont][cont2].a[cont3][cont4]);
									
				}
			}
			printf("Beta = %f ind1 = %d ind2 = %d\n", betar[cont][cont2], cont, cont2);
			if(betar[cont][cont2]<0.001)
			printf("indices de beta que ha dado cero %d %d\n", cont, cont2);
		}	
	 
	 
 }
printf("final de las betas\n");
getchar();
 for(cont=0; cont<NumAng; cont++)
	{
		for(cont2=(factordim-1)*dimx/2; cont2<(factordim+1)*dimx/2; cont2++)
		{
			for(cont3=(factordim-1)*dimy/2; cont3<(factordim+1)*dimy/2; cont3++)
			
			{
				for(cont4=0; cont4<NumPix; cont4++)
				{	
					gamma[cont].a[cont2][cont3]+=(coe[cont][cont4].a[cont2][cont3]); /*Esta bien esta inicializacion¿? */
					printf("%f\n", coe[cont][cont4].a[cont2][cont3]);
					
				}
				
				printf("\n\n");
			//~ printf("gamma=%lf, cont=%d, cont2=%d, cont3=%d\n", *(gamma+cont*dimtotal+cont2*dimy+cont3), cont, cont2, cont3 );
			printf("gamma = %f, ind1 = %d, ind2 = %d, ind3 = %d\n", gamma[cont].a[cont2][cont3], cont, cont2, cont3);
			if(gamma[cont].a[cont2][cont3]<0.0001)
			printf("indices de gamma que ha dado cero %d %d %d y gamma %f\n", cont, cont2, cont3, gamma[cont].a[cont2][cont3]);
			}
			
		}
		
	}
	printf("End of gammas");
getchar();
i=0;
deltaf=1;
printf("comienza iteración\n");
  while((deltaf>0.02) && (i < 200))
 {
 deltaf=0;
 sumpix=0;
	 for(cont=0; cont<NumAng; cont++)
	 {	
		 
		 //~ #pragma omp parallel
		 {
		 
			 //~ #pragma omp for private(cont2, cont3, cont4, cont5, cont6)
		for(cont2=(factordim-1)*dimx/2; cont2<(factordim+1)*dimx/2; cont2++) /*Cambiados limites aqui */
		{//~ printf("cont2=%d, numhilo=%d\n", cont2, omp_get_thread_num());
			for(cont3=(factordim-1)*dimy/2; cont3<(factordim+1)*dimy/2; cont3++)
			{
				for(cont4=0; cont4 < NumPix; cont4++)
				{
					cant=0;
					//~ cont2=con2+(factordim-1)*dimx/2; /* ¿? ¿?¿?¿?¿?¿?¿?*/
					//~ cont3=con3+(factordim-1)*dimy/2; /* ¿?¿?¿?¿*/
	   /*Aqui habia dos for nesteados que he eliminado*/
				
				
				for(cont5=0; cont5<(nonum[cont][cont4].nonnum); cont5++)
				{
					
					//~ cant+=nonum[cont][cont4].as[cont5]*(*(x+(nonum[cont][cont4].pixeles[0][cont5])*dimy+nonum[cont][cont4].pixeles[1][cont5])); /*  sum_{i=0}^{N} a_{ij}*x_{i}  ¿?*/
					cant+=nonum[cont][cont4].as[cont5]*(x[nonum[cont][cont4].pixeles[0][cont5]][nonum[cont][cont4].pixeles[1][cont5]]);
							
				}
				printf("NumAng = %d, cant=%lf, ejex mod=%d, ejey =%d, cont4=%d, coef=%lf, beta = %f x = %f nonum=%f\n", cont, cant, cont2, cont3, cont4, coe[cont][cont2].a[cont3][cont4], betar[cont][cont4],nonum[cont][cont4].as[cont5] ,x[nonum[cont][cont4].pixeles[0][cont5]][nonum[cont][cont4].pixeles[1][cont5]]);
				//~ getchar();
				diferencia+=(b[cont][cont4]-cant)*((coe[cont][cont2].a[cont3][cont4]))/(betar[cont][cont4]); /* Aquí mirar  ---- (b-a*x)*a_{j}/b_j */
				if(diferencia<(-100000))
				{printf("%f", nonum[cont][cont4].as[cont5]*(x[nonum[cont][cont4].pixeles[0][cont5]][nonum[cont][cont4].pixeles[1][cont5]] ));
					getchar();
				}
			} /* Final del bucle que recorre todos los detectores  */
		
	
			x[cont2][cont3]+=lambda*diferencia/(gamma[cont].a[cont2][cont3]); /* update de un pixel */
			deltaf+=diferencia;
			sumpix+=x[cont2][cont3];
			
			
			//~ diferencia=0;
			
			printf("diferencia =%f, pixel nuevo =%f, pixel %d %d, lambda=%f, gamma =%f, beta = %lf, b=%lf\n", lambda*diferencia/(gamma[cont].a[cont2][cont3]), x[cont2][cont3], cont2, cont3, lambda, gamma[cont].a[cont2][cont3], betar[cont][cont4], b[cont][cont4]);
		

	//~ if(diferencia< (-10000))	{printf("eh\n");
		//~ getchar();}		
		
		
		diferencia=0;
		}	/*Final bucle dimy*/
	 //~ getchar();	
	}
}/*Final de parte paralela - aqui he quitado una llave*/

	
 }/*Final de un angulo*/
  deltaf=deltaf/sumpix;
 deltaf=fabs(deltaf);
 /* arg_{min}   */

	 printf("Deltaf=%lf y n de iteraciones =%d\n", deltaf, i);
	 sprintf(filename, "res/Resultados_SART/Resultados_SART_%d.txt", i);
	 
	 finp=fopen(filename, "w+");

for(cont=0; cont<limy; cont++)
{
	
	for(cont2=0; cont2<limx; cont2++)
	{
		
		
		fprintf(finp, "%d %d %f\n",cont, cont2, x[cont][cont2]); 
	}
	
}
//~ *(g+m)+=*(a+m*dimI+i)*(*(p+i)-(*(at+m*dimI+i)*(*(g+ m*dimI+i)))/((*(at+m*dimI+i))*(*(a+m*dimI+i)));
fclose(finp);
	 getchar();
 i+=1;
}
//~ final=clock();
//~ printf("Tardo %lf segundos para %d angulos, %d pixeles y tamaño %d x %d\n", (double)(final-inicio)/CLOCKS_PER_SEC, NumAng, NumPix, dimx, dimy);
//~ finp=fopen("Resultados_SART.txt", "w+");

//~ for(cont=0; cont<dimy; cont++)
//~ {
	
	//~ for(cont2=0; cont2<dimx; cont2++)
	//~ {
		
		
		//~ fprintf(finp, "%f   ", *(x+dimx*cont+cont2)); 
	//~ }
	//~ fprintf(finp, "\n");
//~ }
//~ *(g+m)+=*(a+m*dimI+i)*(*(p+i)-(*(at+m*dimI+i)*(*(g+ m*dimI+i)))/((*(at+m*dimI+i))*(*(a+m*dimI+i)));
//~ fclose(finp);
free(b);
free(x);
/* Liberar memoria*/
FILE *pipeplot=popen("gnuplot -persist", "w");
/*fprintf(pipeplot, "set title /lo que sea/"   */
//~ fprintf(pipeplot, "set pm3d map\n");
//~ fprintf(pipeplot, "set size square 1, 1\n");
fprintf(pipeplot, "plot '/home/pablo/Resultados_C/LFT/res/Resultados_SART/Resultados_SART_%d.txt'w image \n",(i-1));

//~ fprintf(pipeplot, "set term png\n");
//~ fprintf(pipeplot, "set output 'BP_Imagen_Angulos%d_Rayos%d.png'\n", NumAng, ray);
//~ fprintf(pipeplot, "replot\n");
fflush(pipeplot);
fprintf(pipeplot, "quit\n");
pclose(pipeplot);

return 0;
}


void *carga_archivo(int NumAng, int NumPix, float **b)
{
	int cont, cont2, valu;
	
	FILE *finp;
	
 for(cont=0; cont<NumAng; cont++) /*  Cargar b (esto está mal)*/
 {
	 
	 finp=apertura_archivo(cont);
	 
	 for(cont2=0; cont2<NumPix; cont2++)
	 {
	fscanf(finp, "%d", &valu);
	//~ printf("val=%d\n", valu);
	b[cont][cont2]=log(1+valu); /*Esto aquí lo he cambiado*/
	//~ *(b+cont*NumPix+cont2)=valu;
	printf("val=%d, log(val)=%f, cont=%d\n", valu, b[cont][cont2], cont);
	
	 }
	 
	 fclose(finp);
 }
	
	printf("Archivo cargado, presione alguna tecla\n");
	getchar();
}

FILE *apertura_archivo(int contador)
{

char nombrearchivo[25]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "res/angulo_%d.txt", contador);

if(fopen(nombrearchivo, "r")==NULL)
{
	perror("El archivo no existe o no se puede abrir\n");
}

return fopen(nombrearchivo, "r");



}

void saltar_com(FILE *fin) /* En archivo de texto salta a la siguiente línea tras \n*/
{
	char col;
	while(fscanf(fin, "%c", &col))
	{
		if(col=='\n')
		{
			
			break;
		}
		
	}
	
	
}

FILE *apertura_archivo2(int cont1, int cont2)
{
	
	char file_name[40];
	sprintf(file_name, "res/Resultados_SART/Matriz_coef_SART_%d_%d.txt", cont1, cont2);
	
	
	
	return fopen(file_name, "w+");
}

//~ void guarda_matriz_coeficientes(int numang, int numpi, int lix, int liy, void *loq)
//~ {
	
	//~ FILE *fin;
	//~ loq=(struct) coef;
	//~ int cont, cont2, cont3, cont4;
	
	//~ for(cont=0; cont<numang; cont++)
	//~ {
		
		//~ for(cont2=0; cont2<numpi; cont2++)
		//~ {
			//~ fin=apertura_archivo2(cont, cont2);
			
			//~ for(cont3=0; cont3 <lix; cont3++)
			//~ {
				
				//~ for(cont4=0; cont4<limy; cont4++)
				
				//~ {
					
					//~ fprintf(fin, "%d %d %f\n", cont3, cont4, mfeÇ);
					
				//~ }
				
			//~ }
			
		//~ }
		
	//~ }
	
	
	
	
//~ }
