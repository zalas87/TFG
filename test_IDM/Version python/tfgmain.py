#!/usr/bin/env python
# -*- coding: utf-8 -*-


import TrafficModels as TM
import math as mt
import random
import Extra

#Velocidad de entrada
def vin (va, s, p):
	#min(v0, (sqt(ab^2*T^2 - 2*ab*T*va + va^2 - 4*ab*(s0-s)) - ab*T + va)/2)
	#donde ab = 2*sqrt(a*b)
	#va es la velocidad del vehículo que le precede
	return min(p[5],(mt.sqrt(p[3]**2*p[1]**2-2*p[3]*p[1]*va+va**2-4*p[3]*(p[0]-s))-p[3]*p[1]+va)/2);

#Inhomogeneidad:
def Tin (x, p):
	return p[0]*x+p[1]

def getLlegadas(Nmax, dt, lambd):
	arrivaltimes = []

	for i in range(Nmax):
		nextTime = random.expovariate(lambd)
		if nextTime < dt: #adaptación para qeu incluya la llegada, ya que nuestra simulación es discreta
			nextTime = dt
			
		arrivaltimes.append(nextTime) 

	return arrivaltimes


def simulacion(tmax, smax, dt, lambd, params, paramsInho):
	
	numt = int(mt.ceil(tmax/dt + 1)) #Nº máximo de iteracciones
	Nmax = int(mt.ceil(tmax*lambd)) #Nº máximo de posibles llegadas 
	
	#Obtener las llegadas Poissonianas
	arrivaltimes = getLlegadas(Nmax, dt, lambd)
#	arrivaltimes = [0]
#	for i in range(1200):
#		arrivaltimes.append(2.5)
			
	print 'Tiempo promedio de llegada: ' + str(sum(arrivaltimes)/len(arrivaltimes))		
	
            
	llegadas = []	# Para tener un registro de las llegadas a la vía
	stv_gt = [] #lista de listas, contiene toda los estados de los vehículos mientras están circulando por el tramo de la vía, ordenada en tiempo
	A = {} #Hash, estructura de datos de la simulación, Key(id) -> últimos valores del vehículo
	
	#A [id, s, t, v, a, Dx, Dv, s0, T, v0, a, ab, l]
	A[-1] = [-1, 1000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	
	ids = -1
	deletes = 0
	account = 0
	ta = arrivaltimes[account] #Tiempo de llegada del primer vehículo
	complete = False
	la = 0

	for n in range(numt):
		t = dt * n 
		
		#Añadir un nuevo vehículo  
		if ta <= t and (not complete):
			account += 1
			if account <= (len(arrivaltimes) - 1):	#Nuevo tiempo de acceso	
				ta += arrivaltimes[account] 
			else:
				complete = True
				
			if (len(A) == 1):	#Primer Vehículo
				ids += 1
				A[ids] = [ids, 0, t, params[5], 0, A[-1][1], 0, params[0], params[1], params[5], params[2], params[3], params[4]]
				llegadas.append(1)
			elif A[ids][1] > (params[0] +  params[4]): # Resto de vehículos
				ids += 1
				v_temp = vin(A[ids - 1][3],A[ids - 1][1]-params[4], params)
				A[ids] = [ids, 0, t, v_temp, 0, A[ids - 1][1]- params[4], v_temp-A[ids - 1][3], params[0], params[1], params[5], params[2], params[3], params[4]] #mantiene el id para que en la impresión podamos identificarlo
				llegadas.append(1)
			else:
				llegadas.append(0)
				la += 1			
		else:	
			llegadas.append(0)	    


		#Paso de simulacion
		A = TM.idm(A, dt, deletes)

	
		for i in range(len(A) -1 + deletes): #-1 por el vector dummies
			#Inhomogeneidad Timehead (T)
			if (i >= deletes):
				if ((A[i][1] > paramsInho[2]) and (A[i][8] < paramsInho[3])):
					#Para evitar que tras haber pasado sih +lih siga haciendo esta comparacion ya que es constante
					A[i][8] = min(Tin(A[i][1], paramsInho),paramsInho[3])

	    		#Elimina coches de la carretera cuando alcanzan el final 
				if (A[i][1] >= smax):
					del A[i]			
					deletes = i+1
	    
	    	#Guarda en el campo
			if (i >= deletes):
				#Necesario para que no se sobreesbriba el valor en la función idm ya que modificamos el valor al que apunta, por tanto se modifica también el valor de stv_gt
				p = [A[i][0],A[i][1],A[i][2],A[i][3], A[i][4]]
				stv_gt.append(p)

	print 'Nº de llegadas: ' + str(sum(llegadas)) + '\n'

	ids += 1 #ya que python empieza a contar en 0
	return stv_gt, ids



def main():
	
	##DEFINICIÓN DE PÁRAMETROS

	#dt, tiempo máximo, longitud carretera
	dt    = .25         #(s) 
	tmax  = 3000		#(s)
	smax  = 1500        #(m) 
	
	#Parámetros IDM
	s0    = 2            #(m)
	T     = 1.6          #(s)
	a     = 1.5          #(m/s^2)
	b     = 3            #(m/s^2)
	l     = 5            #(m)
	v0    = 120/3.6      #(m/s)
	ab    = 2*mt.sqrt(a*b)  #(m/s^2)

	#Params = [s0 T a ab l v0];
	params = [s0, T, a, ab, l, v0]

	#Parámetros entrada vehículos
	lambd = 3.0         #1 coche por lambda s 

	#Entrada de vehículos:
	lambd = 1/lambd      #(s)
	

	#Parámetros inhomogeneidad
	sih = 500            #Posición (m)
	lih = 200            #Longitud (m)
	Tf  = 2.6            #T final  (s)
	
	#obtenemos los parÃ¡metros que definen la recta de la inhomogeineidad 
	m   = (Tf-T)/lih
	bi   = T-sih*(Tf-T)/lih

	#ParamsInhom = [m bi sih Tf];
	paramsInho = [m, bi, sih, Tf]
	
	stv_gt,ids = simulacion(tmax, smax, dt, lambd, params, paramsInho)

	Extra.compImplementation(stv_gt, ids, tmax, smax, dt, v0)

if __name__ == "__main__":
    main()
