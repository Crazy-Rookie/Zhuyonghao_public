#! /usr/bin/python3
# --* coding UTF-8 *--

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
	0. reciprocal lattice for effective mass
	1. author and copyleft: yonghao_zhu@163.com
	2. date: 2021-05-08
'''

import os
import numpy as np
import math
import linecache

params = {}

params['poscar']  = 'POSCAR'  # standard POSCAR with vasp format

params['high_k']  = [
					 [[0.50, 0.00, 0.00], [0.00, 0.00, 0.00]],
                     [[0.00, 0.00, 0.00], [0.00, 0.50, 0.00]]
                    ]

params['insert_k'] = 40 # in KPOINTS
params['save_']    = True # True or False


def ReciprocalLattice(params):

	if not os.path.isfile(params['poscar']):
		print('No %s! Exiting!' %params['poscar'])
		exit()

	# unit: 0.1nm
	cell_A = [[float(i) for i in linecache.getline(params['poscar'], 3).split()]]
	cell_B = [[float(i) for i in linecache.getline(params['poscar'], 4).split()]]
	cell_C = [[float(i) for i in linecache.getline(params['poscar'], 5).split()]]

	# unit: Bohr; 1 Bohr = 0.5292 * 0.1nm
	cell_A = np.array(cell_A) / 0.5292 
	cell_B = np.array(cell_B) / 0.5292
	cell_C = np.array(cell_C) / 0.5292

	cell_v = np.dot(cell_A,np.cross(cell_B,cell_C).T)

	cell_A_reci = np.abs(np.cross(cell_B,cell_C)/cell_v) * 2 * math.pi
	cell_B_reci = np.abs(np.cross(cell_A,cell_C)/cell_v) * 2 * math.pi
	cell_C_reci = np.abs(np.cross(cell_B,cell_A)/cell_v) * 2 * math.pi

	print('-------------Print-------------------')
	print('reciprocal lattice vectors(in 1/Bohr):')
	cell_reci = np.vstack((cell_A_reci,cell_B_reci,cell_C_reci))
	print(cell_reci)
	print('-------------------------------------')

	high_k = params['high_k']

	for i in range(len(high_k)):

		print('Path %d:' % (i+1) )
		print(high_k[i][0])
		print(high_k[i][1])
		
		k1 = np.array(high_k[i][0])
		k2 = np.array(high_k[i][1])

		k1_reci = np.dot(k1, cell_reci)
		k2_reci = np.dot(k2, cell_reci)

		dis = np.linalg.norm(np.abs(k2_reci - k1_reci), ord = 2)

		dis_step = dis / (params['insert_k'] - 1)

		k1_k2 = np.array([0 + i*dis_step for i in range(params['insert_k'] )])

		if params['save_']:
			file_name = 'path_%s.txt' % (i+1)
			np.savetxt(file_name, k1_k2, fmt='%0.6f')

		print('-------------')


def main(params):

	ReciprocalLattice(params)

if __name__ == '__main__':
	main(params)