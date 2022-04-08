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
import os
import linecache

#read_PROCAR
import read_PROCAR as rp

##################################################################
'''
	deal with HSE06 PROCAR
	output: energy.npy, procar.npy KLABELS, BANDS-HSE06.dat, and PBANDS.dat
'''
EnergyBand_fun = {}
EnergyBand_fun['nbands']        = 48
EnergyBand_fun['effective_kpt'] = 99
EnergyBand_fun['npath']         = 3
EnergyBand_fun['project_bands'] = True # True or False
EnergyBand_fun['pband_roder']   = [[1],[2],[3,4,5],[6]]
EnergyBand_fun['pband_name']    = [['W'],['Mo'],['S'],['Se']]
##################################################################

def EnergyBand(nbands = 48, effective_kpt = 99, npath=3,
	project_bands = False,
	pband_roder = [[1],[2],[3,4,5],[6]],
	pband_name = [['W'],['Mo'],['S'],['Se']]):

	'''
		pband_roder start 1
	'''
	
	nspin = 1
	
	#read PROCAR
	print('----------read PROCAR----------')
	energy, procar = rp.readProcar(nspin=nspin, nbands=nbands)
	
	energy_eff = energy[-1 * effective_kpt:]
	#test#print(energy_eff.shape)

	procar_eff = procar[:, -1 * effective_kpt:, :, :]
	#test#print(procar_eff.shape)
	print('-------------------------------', '\n')

	#read KPOINTS
	print('----------read KPOINTS----------')
	dis = rp.readKPOINTS(npath=npath, effective_kpt=effective_kpt)
	print('--------------------------------')

	#save BAND.dat
	if not os.path.isfile('BANDS-HSE06.dat'):
		with open('BANDS-HSE06.dat', 'w+') as f:
			f.writelines('#K-Path(1/A) Energy-Level(eV)' + '\n')
			f.writelines('# NKPTS & NBANDS:  %s  %s' %(effective_kpt, nbands) + '\n')

			for band in range(nbands):

				title = '# Band-Index    %s' %(band+1)

				f.writelines(title + '\n')

				k_band = energy_eff[:,band]

				if band % 2 == 1:
					for k in range(effective_kpt):
						k_band_ene = '%.6f' % k_band[k]
						k_pos = '%.5f' % dis[k]
						f.writelines('    ' + str(k_pos) + '    ' + str(k_band_ene) + '\n')

				if band % 2 == 0:
					for k in range(effective_kpt):
						k_band_ene = '%.6f' % k_band[effective_kpt - k - 1]
						k_pos = '%.5f' % dis[effective_kpt - k - 1]
						f.writelines('    ' + str(k_pos) + '    ' + str(k_band_ene) + '\n')					

				f.writelines('\n')

	#save KLABELS
	if not os.path.isfile('KLABELS'):
		with open('KLABELS', 'w+') as f:
			f.writelines('K-Label    K-Coordinate in band-structure plots ' + '\n')
			
			dis_ = dis.reshape([npath, -1])
			
			for i in range(npath):
				k_pos = '%.4f' % dis_[i][0]
				f.writelines('1                  %s' % k_pos + '\n')

			k_pos_last = '%.4f' % dis[-1]
			f.writelines('1                  %s' % k_pos_last + '\n')

	#save PBANDS.dat
	
	if not project_bands:
		print('--END--')
		exit()

	print('----------write PBANDS----------')

	for i in range(len(pband_roder)):

		pband_i = pband_roder[i]
		
		name_i = '%s.dat' % pband_name[i][0]

		print('  **write ', name_i, '**')

		pband_i_0 = procar_eff[:, :, :, pband_i[0] - 1]

		if len(pband_i) == 1:
			pband_i = pband_i_0

		if len(pband_i) > 1:

			for x in range(1, len(pband_i)):
				pband_i_0 += procar_eff[:, :, :, pband_i[x] - 1]

			pband_i = pband_i_0
		
		with open(name_i, 'w+') as f:
			f.writelines('#K-Path          Energy    tot' + '\n')
			f.writelines('# NKPTS & NBANDS:  %s  %s' %(effective_kpt, nbands) + '\n')

			for band in range(nbands):

				title = '# Band-Index    %s' %(band+1)

				f.writelines(title + '\n')

				k_band = energy_eff[:,band]

				if band % 2 == 1:

					for k in range(effective_kpt):
						k_band_ene = '%.6f' % k_band[k]
						k_pos = '%.5f' % dis[k]
						k_tot = '%.3f' % pband_i[0, k, band]
						f.writelines('    ' + str(k_pos) + '    ' + str(k_band_ene) + '    ' + str(k_tot) + '\n')

				if band % 2 == 0:
					for k in range(effective_kpt):
						k_band_ene = '%.6f' % k_band[effective_kpt - k - 1]
						k_pos = '%.5f' % dis[effective_kpt - k - 1]
						k_tot = '%.3f' % pband_i[0, effective_kpt - k - 1, band]
						f.writelines('    ' + str(k_pos) + '    ' + str(k_band_ene) + '    ' + str(k_tot) + '\n')

				f.writelines('\n')
	print('--------------------------------')

def main(nbands, effective_kpt, npath, project_bands, pband_roder, pband_name):
	
	EnergyBand(nbands = EnergyBand_fun['nbands'],
		effective_kpt = EnergyBand_fun['effective_kpt'],
		npath         = EnergyBand_fun['npath'],
		project_bands = EnergyBand_fun['project_bands'],
		pband_roder   = EnergyBand_fun['pband_roder'],
		pband_name    = EnergyBand_fun['pband_name'])

if __name__ == '__main__':
	main(
		nbands        = EnergyBand_fun['nbands'],
		effective_kpt = EnergyBand_fun['effective_kpt'],
		npath         = EnergyBand_fun['npath'],
		project_bands = EnergyBand_fun['project_bands'],
		pband_roder   = EnergyBand_fun['pband_roder'],
		pband_name    = EnergyBand_fun['pband_name'])
