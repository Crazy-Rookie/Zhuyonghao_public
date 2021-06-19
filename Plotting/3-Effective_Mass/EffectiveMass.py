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

##################################
import numpy as np
import os
import linecache   
from scipy import optimize
import matplotlib as mpl        
import matplotlib.pyplot as plt 
#################################

import ExtractBand as EB
import FittingBandEdge as FBE
import ReciprocalLattice as RL

################################
# step1: EB
params1              = {}
params1['TorF']      = False # True or False
params1['band_file'] = 'BAND.dat' # BAND.dat or BANDS-HSE06.dat
params1['ikpt']      = 1 # not change
params1['bands']     = [20, 21] # 2 bands, start 1, including
params1['save_']     = True # False or True
params1['show_']     = True  # False or True

if params1['TorF']:
	print('-----Step 1: EB Running-----')
	if os.path.isfile('ExtractBand.txt'):
		print('ExtractBand.txt Exist. Step1 Exiting...\n')
	else:
		print('ExtractBand.py Running...')
		EB.Extract_Band(params1)
		print('ExtractBand.py DONE\n')

#--------------------------------
# step2: RL
params2            = {}
params2['TorF']    = False # True or False
params2['poscar']  = 'POSCAR'  # standard POSCAR with vasp format

params2['high_k']  = [
					 [[0.50, 0.00, 0.00], [0.00, 0.00, 0.00]],
                     [[0.00, 0.00, 0.00], [0.00, 0.50, 0.00]]
                    ]

params2['insert_k'] = 40 # in KPOINTS
params2['save_']    = True # True or False

if params2['TorF']:
	print('-----Step 2: RL Running-----')
	if os.path.isfile('path_1.txt'):
		print('path_1.txt Exist. Step2 Exit...\n')
	else:
		print('ReciprocalLattice.py Running...')
		RL.ReciprocalLattice(params2)
		print('ReciprocalLattice.py DONE\n')

#--------------------------------
# step3: FBE
params3                = {}
params3['TorF']        = False # True or False
params3['ikpt']        = 1 # not change
params3['ExtractBand'] = 'ExtractBand.txt'
params3['insert_k']    = 40 # in KPOINTS
params3['high_k']      = [
					      [[0.50, 0.00, 0.00], [0.00, 0.00, 0.00]],
                          [[0.00, 0.00, 0.00], [0.00, 0.50, 0.00]]
                        ]
params3['bands']       = [20, 21] # 2 bands, start 1, including
params3['ele_filling'] = ['vb', 'cb']
params3['npoints']     = 8
params3['save_']       = False # True
params3['show_']       = False # True or False

if params3['TorF']:
	print('-----Step 3: FBE Running-----')
	if os.path.isfile('band%s_path1.txt' %params3['bands'][0]):
		print('band%s_path1.txt Exist. Step3 Exit...\n' %params3['bands'][0])
	else:
		print('FittingBandEdge.py Running...')
		FBE.FittingBand(params3)
		print('FittingBandEdge.py DONE\n')