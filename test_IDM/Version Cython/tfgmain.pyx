#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cython
cimport cython

import numpy as np
cimport numpy as np

import TrafficModels as TM
import random
import Extra

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


#Velocidad de entrada
cdef vin (DTYPE_t va, DTYPE_t s, np.ndarray[DTYPE_t, ndim=1] p):
	cdef DTYPE_t out	

	out = min(p[5],(np.sqrt(p[3]**2*p[1]**2-2*p[3]*p[1]*va+va**2-4*p[3]*(p[0]-s))-p[3]*p[1]+va)/2);
	return out

cdef Tin (DTYPE_t x, np.ndarray[double, ndim=1] p):
	cdef DTYPE_t out	

	out = p[0]*x+p[1]
	return out


cdef inhomogeneidad(np.ndarray[DTYPE_t, ndim=2] A, np.ndarray[DTYPE_t, ndim=1] p):

	cdef Py_ssize_t rows = A.shape[0] - 1  #Para no contar el vector dummies
	cdef Py_ssize_t row

	if np.size(A,0) >= 2:
		for row in range(0, rows):
			if A[row,1] > p[2] and A[row,8] < p[3]:
				A[row,8] = min(Tin(A[row,1],p),p[3])

	return A


cdef getLlegadas(int Nmax, float  dt, DTYPE_t lambd):
	cdef list arrivaltimes = []
	cdef float nextTime

	for i in range(Nmax):
		nextTime = random.expovariate(lambd)
		if nextTime < dt: #adaptación para qeu incluya la llegada, ya que nuestra simulación es discreta
			nextTime = dt
			
		arrivaltimes.append(nextTime) 

	return arrivaltimes

cdef simulacion(int tmax, int smax, float dt, DTYPE_t lambd, np.ndarray[DTYPE_t, ndim=1] params, np.ndarray[DTYPE_t, ndim=1] paramsInho): 

	cdef Py_ssize_t numt = np.arange(0, tmax, dt).size + 1 
	cdef unsigned int nums = np.ceil(smax/(params[0] + params[4]))
	cdef np.ndarray[DTYPE_t, ndim=1] times

	cdef unsigned int Nmax = int(np.ceil(tmax*lambd)) #Nº máximo de posibles llegadas 

	cdef list arrivaltimes 
#	arrivaltimes = []
#	arrivaltimes.append(0)
	arrivaltimes = getLlegadas(Nmax, dt, lambd)
#	for i in range(1500):
#		arrivaltimes.append(2.5)

	cdef float t
	cdef int complete = 0

	cdef np.ndarray[DTYPE_t, ndim=2] stv_gt = np.empty((numt*lambd*nums,5), dtype=float)
	cdef np.ndarray[long, ndim=1] llegadas = np.zeros(numt, dtype=int)

	print 'Tiempo promedio de llegada: ' + str(sum(arrivaltimes)/len(arrivaltimes))
	 
	#A [id, s, t, v, a, Dx, Dv, s0, T, v0, a, ab, l]
	cdef np.ndarray[DTYPE_t, ndim=2] A = np.array([[-1.0, 1000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]) 
	cdef Py_ssize_t ids, j, k, la, account
	cdef DTYPE_t v_temp, ta

	ids = 0
	j = 0
	k = 0
	la = 0
	account = 0
	ta = arrivaltimes[account]

	times = np.linspace(0,tmax,numt)
	for t in times:
		#Llegadas Poissonianas, parámetro lambda
		if ta <= t and (complete == 0):
			account += 1
			if account <= (len(arrivaltimes) - 1):    #Nuevo tiempo de acceso
				ta += arrivaltimes[account]
			else:
				complete = 1

			if (np.size(A,0) == 1):    #Primer Vehículo
				A = np.vstack((np.array([[ids, 0, t, params[5], 0, A[0,1], 0, params[0], params[1], params[5], params[2], params[3], params[4]]]), A))
				llegadas[k] = 1
			elif (A[0,1] > (params[0] +  params[4])): # Resto de vehículos
				ids += 1
				v_temp = vin(A[0,3],A[0,1]-params[4], params)
				A = np.vstack((np.array([[ids, 0, t, v_temp, 0, A[0,1]- params[4], v_temp-A[0,3],  params[0], params[1], params[5], params[2], params[3], params[4]]]), A))
				llegadas[k] = 1
			else:
				la +=1
		k += 1

		A = TM.idm(A, dt)

		#Inhomogeneidad Timehead (T)
		A = inhomogeneidad(A, paramsInho)

		#Elimina coches de la carretera cuando alcanzan el final 
		A = A[np.logical_or(A[:,1] <= smax, A[:,0] < 0),:]

	    #Guarda en el campo stv_gt todo menos el vector dummies
		stv_gt[j:j+np.size(A,0) - 1,:] = A[0:-1,0:5]
		j = j+np.size(A,0) - 1

	print 'Nº de llegadas: ' + str(sum(llegadas)) + '\n'
	ids += 1 #ya que python empieza a contar en 0
	return stv_gt, ids
		
	
def main():
	#dt, tiempo m�ximo, longitud carretera
	cdef float dt  = .25         #(s) 
	cdef int tmax  = 3000      #(s)
	cdef int smax  = 1500        #(m) 

	#Parametros IDM
	cdef DTYPE_t s0    = 2.0            	#(m)
	cdef DTYPE_t T     = 1.6          	#(s)
	cdef DTYPE_t a     = 1.5          	#(m/s^2)
	cdef DTYPE_t b     = 3.0            	#(m/s^2)
	cdef DTYPE_t ab    = 2*np.sqrt(a*b)  	#(m/s^2)
	cdef DTYPE_t l     = 5.0            	#(m)
	cdef DTYPE_t v0    = 120/3.6     	#(m/s)


	#Params = [s0 T a ab l v0];
	cdef np.ndarray[DTYPE_t, ndim=1] params = np.array([s0, T, a, ab, l, v0])

	#Parámetros entrada vehículos
	cdef DTYPE_t lambd = 3.0          #1 coche por lambda s 

	#Entrada de vehículos:
	lambd = 1/lambd      #(s)
	
	#Parámetros inhomogeneidad
	cdef int sih = 500           	#Posición (m)
	cdef int lih = 200            	#Longitud (m)
	cdef DTYPE_t Tf  = 2.6            #T final  (s)
    
    #obtenemos los parámetros que definen la recta de la inhomogeineidad 
	cdef DTYPE_t m = (Tf-T)/lih
	cdef DTYPE_t bi = T-sih*(Tf-T)/lih	

	cdef Py_ssize_t ids = 0
	
	#ParamsInhom = [m bi sih Tf];
	cdef np.ndarray[DTYPE_t, ndim=1] paramsInho = np.array([m, bi, sih, Tf])
 
	cdef np.ndarray[DTYPE_t, ndim=2] stv_gt

	stv_gt, ids = simulacion(tmax, smax, dt, lambd, params, paramsInho)

	Extra.compImplementation(stv_gt, ids, tmax, smax, dt, v0)





