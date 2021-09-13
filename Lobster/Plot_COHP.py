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
	1. plot coop or cohp with lobster output
	2. DOSCAR.lobster, COHPCAR.lobster
	3. author: yonghao_zhu@163.com
'''

###################################
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import linecache
###################################
params = {}
params['doscar']    = 'DOSCAR.lobster'
params['cohpcar']   = 'COHPCAR.lobster'
params['poscar']    = 'POSCAR'
params['ispin']     = 1 # 1
params['show_']     = True # True or False
params['save_']     = False 
params['linewidth'] = 2
params['y_range']   = [-5, 5, 2] # [start, end, step]
params['x_range']   = [-6, 6, 3] # [start, end, step]
params['pltmode']   = 'cohp' # dos or cohp (-cohp) 
# pltmode = dos
params['show_tdos']   = True
params['ldos_atom']   = [['1', '2', '3-16'], # start 1, type: '1' or '3-16' (including 3 and 16)
			             ['17-32']]
params['ldos_legend'] = ['atoms_1',
			             'atoms_2'] 
params['color']       = ['Red', 'Blue']			      
# pltmode = cohp
params['bond_No.']        = ['No.1', 'No.2', 'No.3'] # start 1
###################################

def ReadPOSCAR(params):
	'''
		read head of POSCAR
		return elements and the numbers
	'''

	input_file = params['poscar']

	elements = linecache.getline(input_file, 6).split()
	numbers = [int(i) for i in linecache.getline(input_file, 7).split()]

	return elements, numbers

def ReadDOSCAR(params = params):

	'''
		Tdos.shape = [energy_point, 2] (energy, Tdos)
		Pdos.shape = [atoms, energy_point, orbitals]
	'''

	input_file = params['doscar']
	
	print('ReadDOSCAR running...')

	# read atoms
	atoms = int(linecache.getline(input_file, 1).strip('\n').split()[0])
	print('atoms =', atoms)

	# read 6th line to get the energy_point
	energy_point = int(linecache.getline(input_file, 6).strip('\n').split()[2])
	print('energy_point =', energy_point)

	# read tdos from 7th to (7+energy_point)th
	Tdos = np.zeros([energy_point,2])
	for i in range(7, 7+energy_point):
		lines = linecache.getline(input_file, i).strip('\n').split()
		energy = float(lines[0]); tdos = float(lines[1])
		Tdos[i-7][0] = energy; Tdos[i-7][1] = tdos
	print('Tdos.shape =', Tdos.shape)

	# read pdos
	states_list = []
	with open(params['doscar']) as f:
		for i in f:
			if 'Z=' in i:
				states_list.append(len(i.split()) - 7)
	#test#print(states_list)

	Pdos = np.zeros([atoms, energy_point, max(states_list)+1])

	for i in range(atoms):
		for j in range(energy_point):
			line_num = (i+1) * (energy_point+1) + 7 + j
			lines = linecache.getline(input_file, line_num).strip('\n').split()
			for s in range(states_list[i]+1):
				Pdos[i, j, s] = float(lines[s])

	print('Pdos.shape =', Pdos.shape, '(Atoms, energy_point, maxmuim orbitals)')

	print('ReadDOSCAR done!')

	return Tdos, Pdos 

def ReadCOHPCAR(params = params):

	'''
		bond_length = ['No.1...', 'No.2...']
		cohp.shape = [energy_point, bonds]
		icohp.shape = [energy_point, bonds]
	'''

	input_file = params['cohpcar']
	
	print('ReadCOHPCAR running...')

	# read 2th line to get the energy_point
	energy_point = int(linecache.getline(input_file, 2).strip('\n').split()[2])
	print('energy_point =', energy_point)

	# read Average
	bond_length = []
	with open(input_file, 'r') as f:
		data = f.readlines()
		for line in data:
			if 'No.' in line:
				bond_length.append(line)
	bonds = len(bond_length)
	print('number of bonds =', bonds)

	# read cohp
	cohp = np.zeros([energy_point, 1+bonds])
	icohp = np.zeros([energy_point, 1+bonds])

	for i in range(bonds+4, bonds+4+energy_point):
		line = linecache.getline(input_file, i).strip('\n').split()
		ene = float(line[0])
		cohp[i-bonds-4, 0] = ene
		icohp[i-bonds-4, 0] = ene
		for j in range(bonds):
			tmp_cohp = line[2*j+1]
			tmp_icohp = line[2*j+2]
			cohp[i-bonds-4, j] = tmp_cohp
			icohp[i-bonds-4, j] = tmp_icohp

	print('cohp.shape =', cohp.shape)
	print('icohp.shape =', icohp.shape)

	print('ReadCOHPCAR done!')

	return bond_length, cohp, icohp

def Plot(params=params):

	Tdos = params['Tdos']
	Pdos = params['Pdos']
	bond_length = params['bond_length']
	cohp = params['cohp']
	icohp = params['icohp']

	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.figure(figsize=(6.5,5))
	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)

	x_ene = Tdos[:, 0]

	if params['pltmode'] == 'dos':
		print('Ploting dos...')

		# plot Tdos
		y_tdos = Tdos[:, 1]
		if params['show_tdos']:
			Tdos_line = plt.plot(x_ene, y_tdos, c='black', zorder=0, lw=params['linewidth'])
			params['ldos_legend'].insert(0, 'TDOS')

		# plot ldos
		ldos_list = []; ldos = []
		for i in params['ldos_atom']:
			tmp = []
			for j in i:
				if '-' in j:
					jj = j.split('-')
					for m in range(int(jj[0]), int(jj[1])+1):
						tmp.append(m)
				else:
					tmp.append(int(j))
			ldos_list.append(tmp)
		
		for i in range(len(ldos_list)):
			pdos_i = np.zeros_like(Pdos)
			
			for j in ldos_list[i]:
				pdos_i[j - 1, :, 1:] = Pdos[j - 1, :, 1:]
				
			ldos_sum_i_projs = np.sum(pdos_i, axis=2)
			ldos_sum_i_atoms = np.sum(ldos_sum_i_projs, axis=0)
			
			ldos.append(ldos_sum_i_atoms)
		
		# plot ldos
		for i in range(len(ldos_list)):
			line_ldos = plt.plot(x_ene, ldos[i], c=params['color'][i], lw=params['linewidth'])
		
		plt.legend(params['ldos_legend'], loc=1, fontsize='large')
		plt.ylabel("DOS ($eV^{-1}$)",size=20)
		plt.xlabel("Energy (eV)",size=20)

	if params['pltmode'] == 'cohp':
		print('Ploting cohp (-cohp)...')

		# cohp
		cohp_sum = np.zeros([cohp.shape[0]])
		for i in params['bond_No.']:
			bond_i = int(i.split('.')[1])
			print(bond_length[bond_i - 1])

			cohp_sum += cohp[:, bond_i]

		cohp_sum = -1*cohp_sum

		cohp_line = plt.plot(x_ene, cohp_sum, c='black', zorder=0, lw=params['linewidth'])
		plt.plot(x_ene, [0 for i in x_ene])
		plt.plot([0 for i in cohp_sum], cohp_sum)

	# set y ticks, range, lable
	y_range = params['y_range']
	plt.ylim(y_range[:2])
	
	length_y = int((max(y_range[:2]) - min(y_range[:2])) / y_range[2])
	plt.yticks([min(y_range[:2]) + i*y_range[2] for i in range(length_y+1)], 
		       [min(y_range[:2]) + i*y_range[2] for i in range(length_y+1)], size=20)

	# set x ticks, range, lable
	x_range = params['x_range']
	plt.xlim(x_range[:2])
	
	length_x = int((max(x_range[:2]) - min(x_range[:2])) / x_range[2])
	plt.xticks([min(x_range[:2]) + i*x_range[2] for i in range(length_x+1)], 
		       [min(x_range[:2]) + i*x_range[2] for i in range(length_x+1)], size=20)

	if params['pltmode'] == 'cohp':
		plt.plot([min(x_range[:2]), max(x_range[:2])], [0, 0], lw=params['linewidth'])
		plt.plot([0, 0], [min(y_range[:2]), max(y_range[:2])], lw=params['linewidth'])

		plt.ylabel("-COHP",size=20)
		plt.xlabel("Energy (eV)",size=20)

	# plot show
	if params['show_']:
		plt.show()

###################################
def main(params):
	if not os.path.isfile(params['doscar']):
		print('No DOSCAR.lobster!')
		print('Exiting...')
		exit()

	if not os.path.isfile(params['cohpcar']):
		print('No COHPCAR.lobster!')
		print('Exiting...')
		exit()

	if params['ispin'] != 1:
		print('please set ispin = 1!')
		print('Exiting...')
		exit()

	if not os.path.isfile(params['poscar']):
		print('No POSCAR!')
		print('Exiting...')
		exit()

	elements, numbers = ReadPOSCAR(params = params)

	print('elements =', elements)
	print('The numbers =', numbers)

	params['elements'] = elements
	params['numbers'] = numbers

	Tdos, Pdos = ReadDOSCAR(params = params)

	print('++++++++++++++++++++++++++++++++')

	bond_length, cohp, icohp = ReadCOHPCAR(params = params)

	params['Tdos'] = Tdos
	params['Pdos'] = Pdos
	params['bond_length'] = bond_length
	params['cohp'] = cohp
	params['icohp'] = icohp

	print('++++++++++++++++++++++++++++++++')

	Plot(params=params)

if __name__ == '__main__':
	main(params)