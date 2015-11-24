#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def idm(A, dt):
    #A [id, s, t, v, a, Dx, Dv, s0, T, v0, a, ab, l]
                
    if np.size(A,0) >= 2:     #Existe un vehÍculo

        #s*(v,Dv) = s0 + max(0, vT + v*Dv/ab) 
        temp_sa = A[0:-1,7] + np.maximum(0, A[0:-1,3]*A[0:-1,8]+A[0:-1,3]*A[0:-1,6]/A[0:-1,11])

        #AceleraciÓn = a * (1 - (v/v0)4 - (s*(v,Dv)/Dx)²)
        A[0:-1,4] = A[0:-1,10]*(1-(A[0:-1,3]/A[0:-1,9])**4-(temp_sa/A[0:-1,5])**2)

        #Integración balística:
        I2 = A[:,3] + A[:,4]*dt >= 0 
        I3 = np.logical_not(I2)       

        #PosiciÓn
		# x = x0 + vt +(1/2)at²
        A[I2,1] = A[I2,1] + A[I2,3]*dt + A[I2,4]*dt**2/2
		# x = x0 - v²/2a
        A[I3,1] = A[I3,1] - (1/2)*A[I3,3]**2/A[I3,4]
        
        #velocidad
		#v = v0 + at
        A[I2,3] = A[I2,3] + A[I2,4]*dt
		#v = 0 
        A[I3,3] = 0  
        
        #Tiempo
        A[0:-1,2] = A[0:-1,2] + dt
        
        #gap y velocidad relativa
        A[0:-1,5] = np.diff(A[:,1]) - A[1:,12]
        A[0:-2,6] = -np.diff(A[0:-1,3])
  
    return A
