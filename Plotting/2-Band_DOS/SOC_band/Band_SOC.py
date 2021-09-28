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

import numpy as np
import os, linecache
#-----------------------------------------
params = {}
# deal with PROCAR within SOC (PBE)
params['poscar'] = 'POSCAR' # not change
params['outcar'] = 'OUTCAR' # not change
params['procar'] = 'PROCAR' # not change
params['BAND']   = 'BAND.dat' # not change, vaspkit format
params['soc_']   = True     # not change
# 1->up; 0->dw (projected z direction)
params['projected_direction'] = 3 # 1->x, 2->y, 3->z
#-----------------------------------------

def readProcar(params):
	'''
	deal with PROCAR within SOC (PBE)
	'''

	proj_dir = params['projected_direction']

	E_f = 0
	with open('OUTCAR') as outcar:
	    for line in outcar.readlines():
		    if 'E-fermi :' in line:
			    lines = line.split()
			    E_f = float(lines[2])

	nkpts = params['nkpts']
	nbands = params['nbands']
	natoms = params['natoms']

	if os.path.isfile('energy.npy'):
		print('Reading the energy.npy...')
		energy = np.load('energy.npy')
	else:
		print('Reading the PROCAR...')
		procar = open('PROCAR','r')
		energy_0 = []
		for i in procar:
			if 'energy' in i:
				energy_0.append(float(i.split()[4])-E_f)
		procar.close()
		energy = np.array(energy_0).reshape(nkpts, nbands)
		np.save('energy.npy',energy)
	print('energy.shape=',energy.shape)

	print('--------------------------------')

	if os.path.isfile('procar.npy'):
		print('Reading the procar.npy...')
		procar = np.loadtxt('procar.npy').reshape(nkpts, nbands, 4, natoms)

	else:
		print('Reading the PROCAR...')
		comments = ['#','band',' k-point','tot','\n','PROCAR','ion',' \n']
		procar_0 = np.loadtxt('PROCAR',comments=comments,dtype=str,delimiter='\n')
		procar_1 = np.array([float(i.split()[-1]) for i in procar_0])
		procar = procar_1.reshape(nkpts, nbands, 4, natoms)
		np.savetxt('procar.npy',procar_1)

	print('PROCAR.shape=',procar.shape)

	print('--------------------------------')

	proj_pro = procar[:, :, proj_dir, :]
	print('Projected spin (directions: %s)' %proj_dir,'[1->x, 2->y, 3->z]')
	print('proj_pro.shape=',proj_pro.shape)
	np.save('proj_pro.npy',proj_pro)

	return proj_pro

def ReWriteBand(params):

	proj_pro = readProcar(params)

	nkpts = params['nkpts']
	nbands = params['nbands']
	natoms = params['natoms']

	print('--------------------------------')

	print('Run ReWriteBand function...')

	reband = np.zeros([nkpts, nbands])

	for ikpt in range(nkpts):
		for iband in range(nbands):
			tot = np.sum(proj_pro[ikpt,iband])
			if tot > 0: # spin up (projected z)
				reband[ikpt,iband] = 1

	reband_list = []
	for iband in range(nbands):
		for ikpt in range(nkpts):
			if iband % 2 == 0: # odd bands, positive sequence 
				reband_list.append(reband[ikpt,iband])
			else: # even bands, negative sequence
				reband_list.append(reband[-1*ikpt-1,iband])

	# rewrite BAND.dat
	up = reband_list # up->1
	dw = [] # dw->1
	for i in up:
		if i == 0:
			dw.append(1)
		else:
			dw.append(0)
	
	data_band = np.loadtxt('BAND.dat')
	
	up = np.array(up); dw = np.array(dw)

	up_dw = np.vstack([up,dw]).T
	
	reband = np.hstack([data_band, up_dw]).reshape(nbands, nkpts, -1)

	with open('ReBand.dat','w+') as f:
		f.writelines('#K-Path(1/A) Energy-Level(eV)'+'\n')
		f.writelines('# NKPTS & NBANDS: %s  %s' % (nkpts, nbands)+'\n')

		for iband in range(nbands):
			f.writelines('# Band-Index    %s' %(iband+1)+'\n')

			for ikpt in range(nkpts):
				for i in reband[iband,ikpt]:
					f.writelines('   '+'%0.5f' % i)
				f.writelines('\n')

			f.writelines('\n')

#-----------------------------------------
def main(params):
	
	#OUTCAR: read E-fermi level
	if not os.path.isfile('OUTCAR'):
		print('No OUTCAR!')
		exit(0)
	#PROCAR: read band energy and projected band 
	if not os.path.isfile('PROCAR'):
		print('No PROCAR')
		eixt(0)

	line_2 = linecache.getline('PROCAR',2).split()
	params['nkpts']  = int(line_2[3])
	params['nbands'] = int(line_2[7])
	params['natoms'] = int(line_2[-1])

	if params['soc_']:
		ReWriteBand(params)
	else:
		print('READ PROCAR with SOC!')

if __name__ == "__main__":
	main(params)
