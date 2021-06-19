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
	1. need BAND.dat (vaspkit -task 211)
	2. author: yonghao_zhu@163.com
	3. date: 2021-01-24
'''

import os                       
import linecache                
import numpy as np             
import matplotlib as mpl        
import matplotlib.pyplot as plt 


params = {}
params['band_file'] = 'BAND.dat' # BAND.dat or BANDS-HSE06.dat
params['ikpt']      = 1 # not change
params['bands']     = [20, 21] # 2 bands, start 1, including
params['save_']     = False # False or True
params['show_']     = True  # False or True

def Extract_Band(params):

	if params['ikpt'] != 1:
		print('ikpt = 1, exitting...')
		exit()

	band_file = params['band_file']
	bands = params['bands']
	nbands = int(linecache.getline(band_file,2).split()[-1])
	save_ = params['save_']
	show_ = params['show_']

	band_matrix = np.loadtxt(band_file).reshape(nbands, -1, 2)

	band_need = np.zeros([len(bands), band_matrix.shape[1], 2])

	for i in range(len(bands)):
		band_need[i, :, :] = band_matrix[bands[i]-1, :, :]
	
	if show_:
		
		plt.rc('font',family='Times New Roman')
		mpl.rcParams['xtick.direction'] = 'in'
		mpl.rcParams['ytick.direction'] = 'in'

		plt.figure(figsize=(6.5,5))
		plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)

		x = band_need[:, :, 0]
		y = band_need[:, :, 1]

		for i in range(len(bands)):

			x_ = x[i]
			y_ = y[i]

			plt.plot(x_, y_, lw = 4, linestyle = ':', marker = 'o')

		plt.ylabel("band energy (eV)",size=20)

		plt.show()

	if save_:

		np.savetxt('ExtractBand.txt', band_need.reshape(-1,2) , fmt = '%.6f')

def main(params):

	Extract_Band(params)

if __name__ == '__main__':
	main(params)