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
	1. allocating nac task within Linux or windows
	2. authour: yonghao_zhu@163.com
'''

import os

###################################################
param = {}
param['time_ini'] = 1 # including, starting 1
param['time_fin'] = 3 # including, starting 1
param['inputs']   = './input_files/' # not change
###################################################

def mkdir(path):

	'''
	make a dir
	'''

	folder = os.path.exists(path)

	if not folder:
		print('mkdir %s...' %path)
		os.makedirs(path)
	
def copy(path_1, path_2):

	'''
	copy path_1 to path_2
	'''

	if not os.path.isfile(path_2):
		with open(path_1, 'r') as f:
			with open(path_2, 'w+') as g:
				for line in f:
					g.writelines(line)

def TaskAllocation(param):

	'''
	time_ini = xxxx
	time_fin = xxxx
	./inputs/(pfiles + KPOINTS + POTCAR + FIRSTCAR + STARTCAR)
	p0001-p9999
	'''
	time_ini, time_fin = param['time_ini'], param['time_fin']
	for time in range(time_ini, time_fin+1):
		print('Time = %s' %time)

		path_tmp = './scf_run/%04d' %time
		path_inputs = param['inputs']
		# mkdir 
		mkdir(path_tmp)
		# INCAR
		if time == time_ini:
			copy(path_1=path_inputs+'FIRSTCAR', path_2=path_tmp+'/INCAR')
		else:
			copy(path_1=path_inputs+'STARTCAR', path_2=path_tmp+'/INCAR')
		# KPOINTS
		copy(path_1=path_inputs+'KPOINTS', path_2=path_tmp+'/KPOINTS')
		# POTCAR
		copy(path_1=path_inputs+'POTCAR', path_2=path_tmp+'/POTCAR')
		# POSCAR
		copy(path_1=path_inputs+'/pfiles/p%04d' %time, path_2=path_tmp+'/POSCAR')

def main(param):

	# check input files
	if not os.path.isfile(param['inputs']+'KPOINTS'):
		print('No %sKPOINTS! Exiting...' % param['inputs'])
		exit()
	elif not os.path.isfile(param['inputs']+'POTCAR'):
		print('No %sPOTCAR! Exiting...' % param['inputs'])
		exit()	
	elif not os.path.isfile(param['inputs']+'FIRSTCAR'):
		print('No %sFIRSTCAR! Exiting...' % param['inputs'])
		exit()	
	elif not os.path.isfile(param['inputs']+'STARTCAR'):
		print('No %sSTARTCAR! Exiting...' % param['inputs'])
		exit()
	elif not os.path.isfile(param['inputs']+'pfiles/p%04d' % param['time_ini']):
		print('No %spfiles/p%04d ! Exiting...' % (param['inputs'],param['time_ini']))
		exit()
	elif not os.path.isfile(param['inputs']+'pfiles/p%04d' % param['time_fin']):
		print('No %spfiles/p%04d ! Exiting...' % (param['inputs'],param['time_fin']))
		exit()
	else:
		print('--------------------------------------------')
		print('Inputs are OK! TaskAllocation Running...')
		TaskAllocation(param)

if __name__ == '__main__':
	main(param)

