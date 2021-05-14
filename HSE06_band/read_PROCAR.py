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
	1. OUTCAR, PROCAR, KPOINTS(hse06), and POSCAR 
	2. no spin 
	3. author: yonghao_zhu@163.com
	4. date: 2021-05-14
'''

import numpy as np
import os, math
import linecache

def readPOSCAR(print_=False):

	atoms_kind = linecache.getline('POSCAR',6).split()
	atoms_num = [int(i) for i in linecache.getline('POSCAR',7).split()]
	#test#print(atoms_kind,atoms_num)

	if print_:
	    print('*******POSCAR*******')
	    print('',atoms_kind,'\n',atoms_num)
	    for i in range(len(atoms_kind)):
	        if i == 0:
	            print('',atoms_kind[i],'-->','1 -',atoms_num[i])
	        if i > 0:
	            print('',atoms_kind[i],'-->',sum(atoms_num[:i])+1,'-',sum(atoms_num[:i])+atoms_num[i])
	    print('********************')

	comments = [linecache.getline('POSCAR',i).rstrip('\n') for i in range(9)]
	del comments[0]
	pos_np = np.loadtxt('POSCAR',comments=comments)[:sum(atoms_num),:]

	return atoms_kind,atoms_num,pos_np
	#test#print(pos_np.shape)	

def readProcar(nspin=1, nbands=792):

	nspin=1 # not change

	#OUTCAR: read E-fermi level
	if not os.path.isfile('OUTCAR'):
		print('No OUTCAR!')
		exit()

	E_f = 0
	with open('OUTCAR') as outcar:
	    for line in outcar.readlines():
		    if 'E-fermi :' in line:
			    lines = line.split()
			    E_f = float(lines[2])
	
	#POSCAR: read the number of atoms
	if not os.path.isfile('POSCAR'):
		print('No POSCAR')
		exit()

	atoms_kind,atoms_num,pos_np = readPOSCAR(print_=False)

	natoms = sum(atoms_num)

	#KPOINTS (hse06): 
	nkpt = int(linecache.getline('KPOINTS', 2))

	#PROCAR: read band energy and projected band 
	if not os.path.isfile('PROCAR'):
		print('No PROCAR')
		eixt()

	if os.path.isfile('energy.npy'):
		energy = np.load('energy.npy')
	else:
		procar = open('PROCAR','r')
		energy_0 = []
		for i in procar:
			if 'energy' in i:
				energy_0.append(float(i.split()[4])-E_f)
		procar.close()
		energy_1 = np.array(energy_0)
		energy_2 = energy_1.reshape(nkpt, nbands)
		np.save('energy.npy',energy_2)
		energy = energy_2
	print('energy.shape=',energy.shape)

	if os.path.isfile('procar.npy'):
		procar = np.loadtxt('procar.npy').reshape(nspin, nkpt, nbands, natoms)
	else:
		comments = ['#','band',' k-point','tot','\n','PROCAR','ion',' \n']
		procar_0 = np.loadtxt('PROCAR',comments=comments,dtype=str,delimiter='\n')
		procar_1 = np.array([float(i.split()[-1]) for i in procar_0])
		procar = procar_1.reshape(nspin, nkpt, nbands, natoms)
		np.savetxt('procar.npy',procar_1)

	print('PROCAR.shape=',procar.shape)

	return energy,procar

def readKPOINTS(npath=3, effective_kpt=99):

	'''
		HSE06 band calculation
	'''

	lattice = []
	for i in range(3,6):
		lines = linecache.getline('POSCAR', i).split()
		cell = [[float(i) for i in lines]]
		lattice.append(cell)
	
	cell_A = np.array(lattice[0])
	cell_B = np.array(lattice[1])
	cell_C = np.array(lattice[2])

	cell_v = np.dot(cell_A,np.cross(cell_B,cell_C).T)
	print('volume of cell= ',round(cell_v[0][0],6),'angstrom^3')

	cell_A_reci = np.abs(np.cross(cell_B,cell_C)/cell_v) * 2*math.pi
	cell_B_reci = np.abs(np.cross(cell_A,cell_C)/cell_v) * 2*math.pi
	cell_C_reci = np.abs(np.cross(cell_B,cell_A)/cell_v) * 2*math.pi
	print('*************************************')
	print('reciprocal lattice vectors(like vasp * 2*pi):')
	cell_reci = np.vstack((cell_A_reci,cell_B_reci,cell_C_reci))
	print(cell_reci)
	print('*************************************')

	with open('KPOINTS', 'r') as f:
		
		k = f.readlines()

		kpoints_list = []
		for i in range(3, len(k)):
			lines = k[i].split()
			if len(lines) == 4:
				tmp = [float(lines[0]), float(lines[1]), float(lines[2])]
				kpoints_list.append(tmp)

		kpoints_list = np.array(kpoints_list)[-1 * effective_kpt:]

	#test#print(kpoints_list)

	Klist = np.dot(kpoints_list, cell_reci).reshape(npath, -1, 3)

	dis = []; dis_0_1 = 0

	for path in range(npath):

		path_0 = Klist[path][0]
		path_1 = Klist[path][-1]

		dis_tmp = []

		for k in Klist[path]:

			dis_ = np.linalg.norm(k - path_0)

			x = dis_0_1 + dis_
			dis.append(x)

		dis_0_1 += np.linalg.norm(path_1 - path_0)

	dis = np.array(dis)

	return dis
