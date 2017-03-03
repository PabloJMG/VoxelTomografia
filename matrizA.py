import math as ma
import numpy as np
import random as rnd
import matplotlib as plt
fdim=2
dim=8
NumPix=8
NumAng=8
limx=limy=dim*fdim
dang=2*ma.pi/NumAng
A=np.zeros((NumPix*NumAng, (fdim*fdim*dim*dim)), dtype=np.float)
fpos=dim/NumPix
def abre_arch(cont1, cont2, lix, liy, *a):
	nom='respy/mat_%d_%d' %(cont1, cont2)
	filen=open(nom, "w")
	for a1 in range(lix):
		for a2 in range(liy):
			
			filen.write(str(a1) + " " + str(a2) + " " + str(a[a1*liy+a2]) + "\n")
	filen.close()		
	return 0
for a1 in range(NumAng):
	beta=dang*a1
	sing=ma.sin(beta)
	cong=ma.cos(beta)
	for a2 in range(NumPix):
		numpi=a1*NumPix+a2
		for a3 in range(dim):
			px=dim/2-a3
			py=dim/2-a2*fpos  #De esto no estoy nada seguro
			px2=px*cong+py*sing
			py2=py*cong-px*sing
			px=px2+fdim*dim/2
			py=py2+fdim*dim/2
			pxt=int(px)
			pyt=int(py)	
			#~ print 'px =', px, 'pxt=', pxt, 'py =', py, 'pyt=', pyt
			#~ raw_input()
			delx=px-pxt
			dely=py-pyt
			
			A[numpi][pxt*fdim*dim+pyt]+=(2-delx-dely)
			#~ print 'a1 =', a1, 'a2 = ', a2, 'px =', pxt, 'py =', pyt, 'valor =', (2-delx-dely), 'A[][]', A[a1*NumPix+a2][pxt*fdim*dim+pyt]
		
			A[numpi][pxt*dim*fdim+pyt+1]+=(1-delx+dely)
			A[numpi][(pxt+1)*dim*fdim+pyt]+=(1-dely+delx)
			A[numpi][(pxt+1)*dim*fdim+pyt+1]+=(delx+dely)
	abre_arch(a1, a2, limx, limy, A[numpi])


print A[2]

#~ def abre_arch(cont1, cont2, lix, liy, *a):
	#~ nom='/respy/mat_%d_%d' %(cont1, cont2)
	#~ filen=open(nom, "w")
	#~ for a1 in range(lix):
		#~ for a2 in range(liy):
			#~ filen.write(cont1, cont2, a[a1*liy+a2])
	#~ filen.close()		
	#~ return 0
	
	
