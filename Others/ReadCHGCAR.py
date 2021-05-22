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
	1. Read CHGCAR, vasp output file, to calculate dipole moment
	2. save ./chgcar/xxxx.npy
	3. author: yonghao_zhu@163.com
	4. date: 2021-05-21
'''

import numpy as np
import linecache

r_f               = {}
r_f['file_g']     = 'CHGCAR'
r_f['valenc_ele'] = [6, 6, 6, 6]
#lattice_z: such as 0 0 25
r_f['direction']  = 'z' # support only z direction, not change
r_f['save_npy']   = False # True or False   


def readCHGCAR(file_g='CHGCAR', direction='z', valenc_ele=[6,6,6,6], save_npy=False):

	'''
		CHGCAR --> 5 data/line
		This only supports z direction!
		the angle between z and x or x and y = 90°
		1 e*Angstrom = 4.80D
	'''

	#atoms
	atoms = [int(i) for i in linecache.getline(file_g,7).split()]
	atoms_num = sum(atoms)

	#lattice_z
	lattice_z = float(linecache.getline(file_g, 5).split()[2])

	#positions_z
	pos = []
	for i in range(atoms_num):
		line = linecache.getline(file_g, i+9).split()
		pos.append(float(line[2]))

	mode = linecache.getline(file_g, 8)

	if 'Direct' in mode:
		pos = [lattice_z * i for i in pos]

	#test#print(pos)

	pos_s = []
	for i in range(len(atoms)):
		if i == 0:
			pos_s.append(pos[ :atoms[i] ])

		if i > 0:
			pos_s.append( pos[ sum(atoms[:i]): sum(atoms[:i])+atoms[i] ] )

	#test#print(pos_s)

	#sum Ziezi
	Ziezi = [ sum(pos_s[i]) * valenc_ele[i] for i in range(len(valenc_ele))]
	#test#print(sum(Ziezi))

	fft = [int(i) for i in linecache.getline(file_g, atoms_num+10).split()]

	total_lines = 0
	if (fft[0] * fft[1] * fft[2]) % 5 ==0:
		total_lines = int((fft[0] * fft[1] * fft[2]) / 5)
	else:
		total_lines = int((fft[0] * fft[1] * fft[2]) / 5) + 1

	if total_lines == 0:
		print('ERROR: Reading fft total_lines')
		exit()

	charge = []

	#read fft data
	for i in range(atoms_num+11, atoms_num+11+total_lines):
		line = linecache.getline(file_g, i)
		lines = line.split()
		for ii in lines:
			charge.append(float(ii)/(fft[0] * fft[1] * fft[2]))

	# dipole moment (e/Angstrom)
	avg = []; x_value = []; dm_ = []

	if direction != 'z':
		print('This only supports z direction!')
		exit()
		
	interval = lattice_z / fft[2]

	for i in range(fft[2]):
		
		start = fft[0] * fft[1] * i
		end = fft[0] * fft[1] * (i+1)
		avg.append(sum(charge[start: end]))

		x_value.append(interval * i)

		dm_i = interval * i * sum(charge[start: end])

		dm_.append(dm_i)

	#test#print(x_value)
	#test#print(dm_); print(sum(dm_))
	
	#save avg
	if save_npy:
		avg_ = np.array(avg).reshape([1,-1])
		np.save('./chgcar/direction_z.npy', avg_)
		
	#save dipole moment
	dm = abs(sum(Ziezi) - sum(dm_))
	with open('DipoleMoment.dat', 'a+') as f:
		f.writelines(str(dm)+' '+'e*Angstrom'+'\n')

	print(abs(sum(Ziezi) - sum(dm_)))


def main():

	readCHGCAR(file_g=r_f['file_g'] , direction=r_f['direction'], valenc_ele=r_f['valenc_ele'],
		       save_npy=r_f['save_npy'])

if __name__ == "__main__":
	main()

'''
CHGCAR format: 
WS2/MoSeS
1.000
 3.182034    0.000002    0.000000
-1.591004    2.755715    0.000000
 0.000000    0.000000   24.202400
 W    Mo   S    Se
 1     1     3     1
Direct
0.331810  0.666004  0.291974
0.998135  0.999190  0.539192
0.665148  0.332672  0.227342
0.665142  0.332670  0.356524
0.331469  0.665857  0.475482
0.331467  0.665857  0.610389

48   48  360
* * * * *
*       *
* * * * *
'''
