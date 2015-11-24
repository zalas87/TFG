#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cython
cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def idm(np.ndarray[DTYPE_t, ndim=2] A, float dt):

	cdef Py_ssize_t rows = A.shape[0] - 1  #Para no contar el vector dummies
	cdef Py_ssize_t row

	cdef DTYPE_t temp_sa

	if np.size(A,0) >= 2: 	#Existe un vehículo
		for row in range(rows-1, -1, -1):

			temp_sa = A[row,7] + max(0, A[row,3]*A[row,8]+A[row,3]*A[row,6]/A[row,11])
			
			#Aceleración = a * (1 - (v/v0)⁴ - (s*(v,Dv)/Ds)²)
			A[row,4] = A[row,10]*(1-(A[row,3]/A[row,9])**4-(temp_sa/A[row,5])**2)


			#Integración balística:
       
			if (A[row,3] + A[row,4]*dt >= 0):
				# x = x0 + vt +(1/2)at²
				A[row,1] = A[row,1] + A[row,3]*dt + A[row,4]*dt**2/2
				# v = v0 + at
				A[row,3] = A[row,3] + A[row,4]*dt
			else:
				# x = x0 - v²/2a
				A[row,1] = A[row,1] - (1/2)*A[row,3]**2/A[row,4]
				# v = 0 
				A[row,3] = 0

			A[row,2] = A[row,2] + dt


			#gap y velocidad relativa
			if (row == rows-1):
				A[row,5] = A[row+1,1] - A[row,1]
				A[row,6] = 0				
			else:
				A[row,5] = A[row+1,1] - A[row,1] - A[row+1,12]
				A[row,6] = A[row,3] - A[row+1,3]
            
			   

	return A	

 
