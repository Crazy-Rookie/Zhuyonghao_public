#! /usr/bin/python3
# -*- conding=UTF-8 -*-

#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /  Zhu  \
#    //  Yong \\
#   //|  Hao  |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__)))  

'''
	1. author: yonghao_zhu@163.com
	2. date: 2021-06-19
'''

##################################
import numpy as np
import os
import linecache   
from scipy import optimize
import matplotlib as mpl        
import matplotlib.pyplot as plt 
##################################

params = {}

params['ikpt']        = 1 # not change
params['ExtractBand'] = 'ExtractBand.txt'
params['insert_k']    = 40 # in KPOINTS
params['high_k']      = [
					      [[0.50, 0.00, 0.00], [0.00, 0.00, 0.00]],
                          [[0.00, 0.00, 0.00], [0.00, 0.50, 0.00]]
                        ]
params['bands']       = [20, 21] # 2 bands, start 1, including
params['ele_filling'] = ['vb', 'cb']
params['npoints']     = 8
params['save_']       = True # True
params['show_']       = False # True or False

#########################################################################################
if params['ikpt'] != 1:
	print('Support ONLY ikpt = 1! Exiting!')
	eixt()

def reorder(ene, insert_k):

	b = np.argsort(ene, axis=0)

	ene_ = np.zeros_like(ene)

	for i in range(params['insert_k']):
		ene_[i] = ene[b[i][0]]

	return ene_


def ExtractPoint(xy, x, npoints, ele_filling):

	print('-----ExtractPoint() function running...-----')

	x_return = []
	y_return = []

	if ele_filling == 'vb':
		
		vbm_v = max(xy[:, 1])
		m_index = int(xy[:, 1].argmax())

	if ele_filling == 'cb':
		
		cbm = min(xy[:, 1])
		m_index = int(xy[:, 1].argmin())

	# right
	if m_index == len(xy[:, 1])-1:

		if ele_filling == 'vb':
			print('The vbm locates at right position!')
		if ele_filling == 'cb':
			print('The cbm locates at right position!')

		y_retrun = xy[:, 1][-1*npoints:]
		x_return = x[-1*npoints:]

	# left
	if m_index == 0:

		if ele_filling == 'vb':
			print('The vbm locates at left position!')
		if ele_filling == 'cb':
			print('The cbm locates at left position!')

		y_retrun = xy[:, 1][:npoints]
		x_return = x[:npoints]

	# middle
	if m_index > 0 and m_index < len(xy[:, 1])-1:

		if ele_filling == 'vb':
			print('The vbm locates at middle position!')
		if ele_filling == 'cb':
			print('The cbm locates at middle position!')

		if npoints % 2 == 0:
			npoints_ = int(npoints/2)
		else:
			npoints_ = int(npoints) + 1
		y_retrun = xy[:, 1][m_index-npoints_:m_index+npoints_]
		x_return = x[m_index-npoints_:m_index+npoints_]


	print('x_return:', x_return)
	print('y_retrun:', y_retrun)

	print('--------------------DONE--------------------')

	return x_return, y_retrun


def QuadraticFunction(x, a, b, c):
	return a* x**2 + b*x + c 

def R2(a, b, c, x, y_real):
	'''
	y_ = sum(yi)/n
	SStot = sum((yi-y_)^2)
	SSreg = sum((fi-y_)^2)
	SSres = sum((yi-fi)^2)
	R^2 = 1- SSres/SStot = SSreg/SStot
	'''

	y_fitting = [ a*i**2 + b*i + c for i in x]

	y_ = sum(y_real) / len(y_real)

	SSreg = sum([ (y_fitting[i] - y_)**2 for i in range(len(x)) ])

	SStot = sum([ (y_real[i] - y_)**2 for i in range(len(x)) ])

	print('R^2 = %s' %(SSreg/SStot))


def FittingBand(params):

	npath = len(params['high_k'])
	nbands = len(params['bands'])

	ene = np.loadtxt(params['ExtractBand']).reshape(nbands, npath, params['insert_k'], 2)
	#test#print(ene)

	if params['save_'] and not os.path.isfile('band%s_path%s.dat' % (params['bands'][0], 1)):
		print('-----Save bandVBM_path1.dat-----')
		
		for band in range(1, nbands+1):
			if band % 2 == 0:
				for path in range(1, npath+1):
					file_name = 'band%s_path%s.dat' % (params['bands'][band-1], path)
					np.savetxt(file_name, reorder(ene[band-1, path-1], insert_k=params['insert_k']), fmt='%.4f')

			if band % 2 == 1:
				for path in range(1, npath+1):
					file_name = 'band%s_path%s.dat' % (params['bands'][band-1], npath-path+1)
					np.savetxt(file_name, reorder(ene[band-1, path-1], insert_k=params['insert_k']), fmt='%.4f')				


		print('--------------DONE--------------')


	# Fitting
	print('---------------Fitting Band Edge---------------\n')
	for band in range(1, nbands+1):

		ele_filling = params['ele_filling'][band-1]

		for path in range(1, npath+1):
			
			file_ene = 'band%s_path%s.dat' % (params['bands'][band-1], path)
			file_path = 'path_%s.txt' % path

			# ev --> Hartree
			xy = np.loadtxt(file_ene)*0.036749
			x = np.loadtxt(file_path)

			if ele_filling == 'vb':
				x,y = ExtractPoint(xy=xy, x=x, npoints=6, ele_filling='vb')

			if ele_filling == 'cb':
				x,y = ExtractPoint(xy=xy, x=x, npoints=6, ele_filling='cb')

			print('energy file: %s' %file_ene)
			print('path file: %s' %file_path)

			#fitting
			popt, pcov = optimize.curve_fit(QuadraticFunction, x, y)
			
			a = popt[0]
			b = popt[1]
			c = popt[2]
			
			R2(a=a, b=b, c=c, x=x, y_real=y)

			print('Fitting a=%s' %a)
			print('Effective Mass(m0)=%s' % (0.5/a))

			if params['show_']:
				plt.rc('font',family='Times New Roman')
				mpl.rcParams['xtick.direction'] = 'in'
				mpl.rcParams['ytick.direction'] = 'in'
				plt.figure(figsize=(6.5,5))
				plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)
				plt.plot(x, y, lw = 2, linestyle = ':', marker = 'o')
				plt.show()

			print('=============================================================================\n')
			
	print('------------DONE------------')


def mian(params):
	FittingBand(params)

if __name__ == '__main__':
	main(params)
