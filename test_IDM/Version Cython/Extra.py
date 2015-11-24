#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pylab
import xlrd 


def compImplementation(stv_gt, ids, tmax, smax, dt, v0):

	#Dimensión inicial de la implementación python
	treco_py= np.zeros((ids)) 
	
	#Leer excel
	book = xlrd.open_workbook("../Matlab/a9.xls") 
	sh = book.sheet_by_index(0)
	treco_ma = np.zeros((sh.nrows))

	srecmax = v0*dt
	dimension = 0
	n = 0

	#Convierte stv_gt en una array numpy, en el caso de la implementación python
	if type(stv_gt).__module__ != np.__name__:
		stv_gt = np.asarray(stv_gt)
	else: 
		stv_gt = stv_gt[stv_gt[:,1] > 0.0]

	for i in range(ids):
		I = np.nonzero(stv_gt[:,0] == i)


		#devuelve una trupla ...
		entrada =  I[0][0]
		salida = I[0][-1]

		if (stv_gt[salida,1] >= (smax - srecmax)):
			duracion = stv_gt[salida,2] - (stv_gt[entrada,2] - dt)
			treco_py[i] = duracion 

	
    #reajustar dimensión -> Eliminación de 0's
	treco_py = treco_py[treco_py[:] > 0.0]
	
	for row in range(sh.nrows):
		cell = sh.cell(row,0)
		treco_ma[row] = cell.value

	tmean_ma = np.mean(treco_ma) 
	tstd_ma = np.std(treco_ma)
	tmean_py = np.mean(treco_py)
	tstd_py = np.std(treco_py)
	print 'Media Matlab: ' + str(tmean_ma)
	print 'Desviación típica Matlab: ' + str(tstd_ma)
	print 'Media implementación propia: ' + str(tmean_py)
	print 'Desviación típica implementación propia: ' + str(tstd_py) + '\n'

	#Nos quedamos con la simulación que tenga menos vehículos
	if np.size(treco_py) >= np.size(treco_ma):
		dimension = np.size(treco_ma)
	else:
		dimension = np.size(treco_py)
		
	#Promediamos cada 5 vehículos
	mean_py = np.zeros((dimension/5))
	mean_ma = np.zeros((dimension/5))
	std_py = np.zeros((dimension/5))
	std_ma = np.zeros((dimension/5))
	
	ite = np.linspace(5,dimension,dimension/5)
	
	#Contruimos el vector que contiene el valor de tiempo promedio en recorrer el tramo de la vía
	for i in ite:
		mean_py[n] = np.mean(treco_py[(i-5):i]) 
		mean_ma[n] = np.mean(treco_ma[(i-5):i])
		std_py[n] = np.std(treco_py[(i-5):i])
		std_ma[n] = np.std(treco_ma[(i-5):i])
		n += 1

	ideal = np.linspace(min(mean_py) - 5,max(mean_py) + 5,50)
	

	pylab.figure(0)
	pylab.subplot(1,2,1)
	pylab.plot(mean_py,mean_ma,'b.', label = 'Real');
	pylab.title(r"$Matlab \ VS \ Cython:\ T'=2.6,\ \lambda=3.0$", fontsize=12)
	pylab.plot(ideal, ideal, 'k--', label = 'Ideal');
	pylab.legend(loc = 'upper left')
	pylab.xlabel(r'$\bar t_{cython}$');
	pylab.ylabel(r'$\bar t_{matlab}$');

	pylab.subplot(1,2,2)	
	pylab.plot(std_py,std_ma,'b.');
	pylab.title(r"$Matlab \ VS \ Cython:\ T'=2.6,\ \lambda=3.0$", fontsize=12)
	pylab.xlabel(r'$\sigma_{cython}$');
	pylab.ylabel(r'$\sigma_{matlab}$');
	
	pylab.figure(1)
	pylab.scatter(stv_gt[:,2],stv_gt[:,1],30,stv_gt[:,3],marker='.',edgecolors='none');
	pylab.title(r"$Cython:\ T'=2.6,\ \lambda=3.0$", fontsize=12)
	pylab.axis([0, tmax, 0, smax]);
	pylab.xlabel('Tiempo (s)');
	pylab.ylabel('Espacio (m)');
	pylab.colorbar()

	pylab.show()	
	

	
	return 0

