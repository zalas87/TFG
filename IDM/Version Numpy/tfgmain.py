#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import TrafficModels as TM
import random

#Velocidad de entrada
def vin (va, s, p):
    #min(v0, (sqt(ab^2*T^2 - 2*ab*T*va + va^2 - 4*ab*(s0-s)) - ab*T + va)/2)
    return min(p[5],(np.sqrt(p[3]**2*p[1]**2-2*p[3]*p[1]*va+va**2-4*p[3]*(p[0]-s))-p[3]*p[1]+va)/2);

#Inhomogeneidad:
def Tin (x, p):
    return p[0]*x+p[1]

def getLlegadas(Nmax, dt, lambd):
	arrivaltimes = []

	for i in range(Nmax):
		nextTime = random.expovariate(lambd)
		if nextTime <= dt: #adaptación para qeu incluya la llegada, ya que nuestra simulación es discreta
			nextTime = dt
			
		arrivaltimes.append(nextTime) 

	return arrivaltimes

def simulacion(tmax, smax, dt, lambd, params, paramsInho):
    ids = 0

    numt = np.arange(0, tmax, dt).size + 1 
    nums = np.ceil(smax/(params[0] + params[4]))
    
    Nmax = int(np.ceil(tmax*lambd)) #Nº máximo de posibles llegadas 
    
    arrivaltimes = getLlegadas(Nmax, dt, lambd)
                
    llegadas = np.zeros((numt))
    
    stv_gt = np.zeros((numt*nums*lambd,5))

    A = np.array([[-1, 1000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    j = 0
    k = 0
    
    account = 0
    ta = arrivaltimes[account]
    complete = False

    times = np.linspace(0,tmax,numt)

    for t in times:
        #Llegadas Poissonianas, parámetro lambda             
        if ta <= t and (not complete):
            account += 1
            if account <= (len(arrivaltimes) - 1):    #Nuevo tiempo de acceso
                ta += arrivaltimes[account] 
            else:
                complete = True
                
            if (np.size(A,0) == 1):    #Primer Vehículo
                A = np.vstack((np.array([[ids, 0, t, params[5], 0, A[0,1], 0, params[0], params[1], params[5], params[2], params[3], params[4]]]), A))
                llegadas[k] = 1
            elif (A[0,1] > (params[0] +  params[4])): # Resto de vehículos
                ids += 1
                v_temp = vin(A[0,3],A[0,1]-params[4], params)
                A = np.vstack((np.array([[ids, 0, t, v_temp, 0, A[0,1]- params[4], v_temp-A[0,3],  params[0], params[1], params[5], params[2], params[3], params[4]]]), A))
                llegadas[k] = 1
                
        k += 1
        
        #Paso de simulación
        A = TM.idm(A, dt)
 #       print A
        #Inhomogeneidad Timehead (T)
        I1 = A[:,1] > paramsInho[2];
        A[I1,8] = np.minimum(Tin(A[I1,1], paramsInho),paramsInho[3])
        
        #Elimina coches de la carretera cuando alcanzan el final 
        A = A[np.logical_or(A[:,1] <= smax, A[:,0] < 0),:]
        
        #Guarda en el campo
        stv_gt[j:j+np.size(A,0)-1,:] = A[0:-1,0:5]
        j = j+np.size(A,0) - 1

    return stv_gt



def main():
    
    ##DEFINICIÓN DE PARÁMETROS

    #dt, tiempo máximo, longitud carretera
    dt    = .25         #(s) 
    tmax  = 3000       #(s)
    smax  = 1500        #(m) 
    
    #Parámetros IDM
    s0    = 2            #(m)
    T     = 1.6          #(s)
    a     = 1.5          #(m/s^2)
    b     = 3            #(m/s^2)
    l     = 5            #(m)
    v0    = 120/3.6      #(m/s)
    ab    = 2*np.sqrt(a*b)  #(m/s^2)
    
    #Params = [s0 T a ab l v0];
    params = np.array([s0, T, a, ab, l, v0])

    #Parámetros entrada vehículos
    lambd = 2.7         #1 coche por lambda s 

    #Entrada de vehículos:
    lambd = 1/lambd      #(s)

    #Parámetros inhomogeneidad
    sih = 500            #Posición (m)
    lih = 200            #Longitud (m)
    Tf  = 2.6            #T final  (s)
    
    #obtenemos los parámetros que definen la recta de la inhomogeineidad 
    m   = (Tf-T)/lih
    bi   = T-sih*(Tf-T)/lih


    #ParamsInhom = [m bi sih Tf];
    paramsInho = np.array([m, bi, sih, Tf])
    
    stv_gt = simulacion(tmax, smax, dt, lambd, params, paramsInho)

#    np.savetxt('stv_gt.txt', stv_gt, delimiter=' | ',  fmt='%1.4e')

if __name__ == "__main__":
    main()





