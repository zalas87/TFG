#!/usr/bin/env python
# -*- coding: utf-8 -*-

def idm(A, dt, deletes):
    #A [id, s, t, v, a, Dx, Dv, s0, T, v0, a, ab, l]

    if len(A) >= 2:     #Existe un veh�culo
        for i in range(len(A) -1):
            j = i + deletes
        
            #s*(v,Dv) = s0 + max(0, vT + v*Dv/ab)
            temp_sa = A[j][7] + max(0, A[j][3]*A[j][8]+A[j][3]*A[j][6]/A[j][11])
            #Aceleraci�n = a * (1 - (v/v0)^4 - (s*(v,Dv)/Dx)�)
            A[j][4] = A[j][10]*(1-(A[j][3]/A[j][9])**4-(temp_sa/A[j][5])**2)

            #Integración balística:
       
            if (A[j][3] + A[j][4]*dt >= 0):
				# x = x0 + vt +(1/2)at²
				A[j][1] = A[j][1] + A[j][3]*dt + A[j][4]*dt**2/2
				# v = v0 + at
				A[j][3] = A[j][3] + A[j][4]*dt
            else:
				# x = x0 - v²/2a
				A[j][1] = A[j][1] - (1/2)*A[j][3]**2/A[j][4]
				# v = 0 
				A[j][3] = 0


            #Tiempo
            A[j][2] = A[j][2] + dt
            
             #gap y velocidad relativa
            if i == 0: #Si es el leader, ya que el último elemento es del hash es -1 (dummies)
                A[j][5] = A[-1][1] - A[j][1]
                A[j][6] = 0
            else:
                A[j][5] = A[j-1][1] - A[j][1] - A[j-1][12]
                A[j][6] = A[j][3] - A[j-1][3]

    return A
