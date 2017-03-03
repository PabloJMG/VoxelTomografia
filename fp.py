#!/usr/bin/env python
 #-*- coding: utf-8 -*-

import numpy as np
import math as ma
import matplotlib as plt
import random as rnd

#~ import prop
part = {'NRay': 10000,'NumPix':20, 'dimx':100 ,'dimy':100,'meshsize':1.0,'fdim':3,'d_src_mesh':0.0002, 'srcsize':0.005, 'dx':0.1,'dy':0.1, 'detsize':1.0}

print'Numero de angulos  y presiona start: '
NumAng=input()
deltang=2*ma.pi/NumAng
part["dang"]=deltang
limx=part["dimx"]*part["fdim"]
limy=part["dimy"]*part["fdim"]
part['limx']=limx
part['limy']=limx
part['amax']=part["meshsize"]/(2*(part["meshsize"]+part["d_src_mesh"])) #Cambiado porque si
deltx=part["dimx"]*part["fdim"]/10
delty=part["dimy"]*part["fdim"]/10
rad=part["dimx"]*part["fdim"]/10
objeto=np.zeros([limx,limy])
part['Niter']=int(part["dimx"]/part['dx'] )#Aqui cambie
part['deltaP']=part['detsize']/part['NumPix']
#func(**data) interpreta argumento como dict
#func(*data ) " 		"		"	como tupla
print part['amax'], part['meshsize']

param=part
mutot=0
#~ for cont in range(limx):	#Inicializa objeto (cuadrado)
	#~ for cont2 in range(limy):
		#~ if(((math.abs(cont-limx/2) < deltx)) and ((math.abs(cont-limy/2) < delty))) :
		#~ objeto[cont,cont2]=5
		
		#Inicializa objeto (circulo)
for cont1 in range(limx):
	for cont2 in range(limy):
		if(((cont1-limx/2)**2 +(cont2- limy/2)**2)<rad**2):
			objeto[cont1,cont2]=500000
		

for con in range(NumAng):

	
	arpi=np.zeros(part["NumPix"])
	dc2=param["meshsize"]/param["Niter"] #Cambié la definición del dx. Parece que funciona mejor
	sinang=ma.sin(con*part["dang"])
	cosang=ma.cos(con*part["dang"])	
	hit=0
	for cont1 in range(part["NRay"]):
		amax=param["amax"]
		ang=-amax+2*amax*rnd.random()
		pos=-param["srcsize"]+2*param["srcsize"]*rnd.random()
		pos+=ang*param["d_src_mesh"]
		mutot=0
		for cont2 in range(param["Niter"]): #Ni idea de si existe Niter ahora
			pos+=ang*dc2
			pix=-param["dimx"]/2+cont2*dc2 #Cambio aqui
			piy=(pos/param["dy"])
			pix2 = pix*cosang+piy*sinang
			piy2 = piy*cosang-pix*sinang
			pix=int(round(pix2+param["dimx"]*param["fdim"]/2))
			piy=int(round(piy2+param["dimy"]*param["fdim"]/2))

			mu=objeto[pix, piy]
			#~ print 'mu =', mu, 'rayo =', cont1, 'angulo inicial =', ang
			mutot+=mu/param["Niter"]
	
		#~ print 'Mu total =', mutot, 'y exp =',ma.exp(-mutot)	
	
		if(rnd.random()<ma.exp(-mutot*param["dimx"])):
			npix=int(round((pos+param["detsize"]/2)/param["deltaP"])) #npix esta mal
			#~ print pos, (pos+param["detsize"]/2)/param["deltaP"], round((pos+param["detsize"]/2)/param["deltaP"]) , npix
			if(npix>=0 and npix<param["NumPix"]):
				arpi[npix]+=1
				hit+=1
				#~ print 'hola2'
	
	
	filename = 'respy/angulo_%d.txt' %(con)
	fout=open(filename, "w")
	fout.write('\n'.join('%d' % x for x in arpi))
	print arpi, 'Rayos que llegan ',  hit
	#~ raw_input()
	fout.close()




