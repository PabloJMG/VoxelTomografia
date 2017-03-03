
#!/usr/bin/env python
 #-*- coding: utf-8 -*-

import math as ma
import numpy as np
import random as rnd
from scipy.sparse import sparse_cscmatrix
NumAng=10
NumPix=200
dim=64
dimtotal=dim**2
fdim=2
fpos=dim/NumPix
dang=2*ma.pi/NumAng
#Carga resultados de la propaga en b
for con1 in range(NumAng):
	finp=abrir_archivo(con1)
	
	c=finp.readlines()
	
	c=(x.strip() for x in c)
	b.append(c)
	
	#~ print b

for x in range(len(b)):
	b[x]=tuple(b[x])

#~ [i[0] for i in b]	

b=list(sum(b,()))
#~ b=np.zeros[NumAng*NumPix]
x=np.zeros[(fdim*dim)**2]
A=np.zeros[NumAng*NumPix, (fdim*dim)**2] 
#Carga de los archivos a b
#Creacion de la matriz A
	
for a1 in range(NumAng):
	beta=dang*a1
	sing=ma.sin(beta)
	cong=ma.cos(beta)
	for a2 in range(Numpix):
		for a3 in range(dim):
			px=dim/2-a3
			py=dim/2-a2*fpos
			px2=px*cong+py*sing
			py2=py*cong-px*sing
			px=px2+fdim*dim/2
			py=py2+fdim*dim/2
			pxt=int(round(px))
			pyt=int(round(py))	
			delx=px-pxt
			dely=py-pyt
			A[a1*NumPix+a2][px*dim+py]+=(2-delx-dely)
			A[a1*NumPix+a2][px*dim+py+1]+=(1-delx+dely)
			A[a1*NumPix+a2][(px+1)*dim+py]+=(1-dely+delx)
			A[a1*NumPix+a2][(px+1)*dim+py+1]+=(delx+dely)
#Crear matriz sparse
conta=0
sizA = A.shape #En realidad, se podria definir manualmente, nunca cambia su taman en teoria.
C=sparse_cscmatrix(1.0/sum(A), range(sizA[1]),range(sizA[1]))
R=sparse_cscmatrix(1.0/sum(A), range(sizA[0]), range(sizA[0]))
while((inc>0.1)and conta<100):	
	inc=C*AT*R*(b-A*x)
	x+=inc
	
	
def abrir_archivo(con):
	filename = 'respy/angulo_%d.txt' %(con)
	return open(filename, "r")
