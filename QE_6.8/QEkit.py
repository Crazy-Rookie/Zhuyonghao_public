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
	http://www.quantum-espresso.org/
	QEkit (QE6.8)
'''

print('Supportting for QE6.8 version.')

######################################
import sys, os, math
import linecache
import numpy as np
######################################

VASPtoXYZ_fun = {}
'''
    input file:  POSCAR.vasp (vasp format)
            xxx-->character string
            1.0
               1 0 0
               0 1 0
               0 0 1
            m n
            1 1 
            D or C
            0.0 0.0 0.0 
            0.5 0.5 0.5
    output file: POSCAR.xyz (xyz format)
'''
VASPtoXYZ_fun['TorF'] = False # True or False

#-------------------------------------

GetFinPOS_fun             = {}
GetFinPOS_fun['TorF']     = True # True or False
GetFinPOS_fun['log_file'] = 'cell_opt.log'

#-------------------------------------

Getband_fun = {}
Getband_fun['TorF']     = False # True or False
Getband_fun['log_file'] = 'band.log' 
Getband_fun['fermi']    = -8.324160481667850e-2 # au, grep fermi pw_band.xml
Getband_fun['K_lable']  = ['G', 'M', 'K', 'G']
Getband_fun['K_coord']  = [[0.00, 0.00, 0.00],
						   [0.50, 0.00, 0.00],
						   [0.3333, 0.3333, 0.00],
						   [0.00, 0.00, 0.00],]

#-------------------------------------

Getwfn_fun = {}
# get POSCAR.vasp, copy the infor into wfn files
Getwfn_fun['TorF']      = False # True or False 
Getwfn_fun['inputs']    = 'band.in' # band.in or pw.scf (atomic_positions: crystal)
Getwfn_fun['atom_kind'] = 'Mo S'
Getwfn_fun['atom_numb'] = '1 2'

######################################
def VASPtoXYZ(file_VASP='POSCAR.vasp', file_XYZ='qe.geo'):

	#read atoms kinds and number
	atoms_kind = linecache.getline(file_VASP,6).split()
	atoms_num = [int(i) for i in linecache.getline(file_VASP,7).split()]

	#write xyz file
	xyz_file = open(file_XYZ,'w+')
	xyz_file.writelines(str(sum(atoms_num))+'\n')
	xyz_file.writelines('POSCAR.xyz'+'\n')

	#read positions
	comments = [linecache.getline(file_VASP,i).rstrip('\n') for i in range(9)]
	del comments[0]
	print(comments)
	pos_np = np.loadtxt(file_VASP, comments=comments)
	print('POSCAR shape:',pos_np.shape)

	if 'D' in linecache.getline(file_VASP,8): #Direct
	    print('Direct mode is OK!')

	if 'C' in linecache.getline(file_VASP,8): #Cartesian
	    print('Not support Cartesian. Exiting...')
	    exit()

	#write pos_np file
	print('********************')
	for i in range(len(atoms_kind)):
	    if i == 0:
	        print('',atoms_kind[i],'-->','1 -',atoms_num[i])
	        for ii in range(0,atoms_num[i]):
	            xyz_file.writelines(atoms_kind[i]+'    '+
	                '%.9f' % pos_np[ii][0]+'  '+
	                '%.9f' % pos_np[ii][1]+'  '+
	                '%.9f' % pos_np[ii][2]+'  '+'\n')
	    if i > 0:
	        print('',atoms_kind[i],'-->',sum(atoms_num[:i])+1,'-',sum(atoms_num[:i])+atoms_num[i])
	        for ii in range(sum(atoms_num[:i]),sum(atoms_num[:i])+atoms_num[i]):
	            xyz_file.writelines(atoms_kind[i]+'    '+
	                '%.9f' % pos_np[ii][0]+'  '+
	                '%.9f' % pos_np[ii][1]+'  '+
	                '%.9f' % pos_np[ii][2]+'\n')
	print('********************')

if VASPtoXYZ_fun['TorF']:
	print('**********Running VASPtoXYZ!!!**********')
	if os.path.isfile('POSCAR.vasp'):
		VASPtoXYZ(file_VASP='POSCAR.vasp', file_XYZ='qe.geo')
	else:
		print('No POSCAR.vasp!!!')

######################################
def GetFinPOS(logfile, outfile):
	'''
	get CONTCAR from logfile.
	'''
	num_start = 0; num_end = 0; contcar = []

	with open(logfile, 'r') as f:

		for num,line in enumerate(f):
			if 'Begin final coordinates' in line:
				num_start = num
				print(num_start)
			if 'End final coordinates' in line:
				num_end = num
				print(num_end)

		if num_end == num_start or num_end < num_start:
			print('Error num_start or num_end! Exitting...')
			exit()

	with open(outfile, 'w+') as f:
		for i in range(num_start+1, num_end+1, 1):
			line = linecache.getline(logfile, i)
			f.writelines(line)
		f.writelines('=============================================================='+'\n')

		position = []
		# positions writing
		for i in range(num_start+1+10, num_end+1, 1):
			pos = linecache.getline(logfile, i)
			position.append(pos)
			f.writelines(pos[:4])		

		f.writelines('\n'+'==============================================================')	
		f.writelines('\n'+'CONTCAR'+'\n'+'1'+'\n')

		# lattice writing
		for i in range(num_start+1+5, num_start+1+8, 1):
			latt = linecache.getline(logfile, i)
			f.writelines(latt)

		f.writelines('\n'+'\n'+'Direct'+'\n')
		for i in position:
			f.writelines(i[10:])

if GetFinPOS_fun['TorF']:
	print('**********Running GetFunPos!!!**********')
	if os.path.isfile(GetFinPOS_fun['log_file']):
		GetFinPOS(logfile=GetFinPOS_fun['log_file'], outfile='qe_contcar.geo')
	else:
		print('No %s! Exitting...' % GetFinPOS_fun['log_file'])

######################################
def Getband(logfile, fermi, K_lable, K_coord):

	k_coordinates_c = []
	k_energy = []; k_energy_lines = []
	lattice_r = []

	au = 2.72113838565563E+01

	with open(logfile, 'r') as f:
		data = f.readlines()
		for line in data:

			# read lattice matrix
			if 'reciprocal axes:' in line:
				tmp = data.index(line)
				b1 = [float(data[tmp+1].split()[i+3]) for i in range(3)]
				b2 = [float(data[tmp+2].split()[i+3]) for i in range(3)]
				b3 = [float(data[tmp+3].split()[i+3]) for i in range(3)]
				lattice_r.append(b1); lattice_r.append(b2); lattice_r.append(b3)

			# read k point coordinates
			if 'k(   ' in line:
				lines = line.split()
				k_coordinates_c.append([float(lines[4]), float(lines[5]), float(lines[6][:-2])])

			# read k point energy
			if 'bands (ev):' in line:
				k_energy_lines.append(data.index(line))

			if 'Writing output data file' in line:
				k_energy_lines.append(data.index(line))

	k_energy_lines[-2] = k_energy_lines[-1] - (k_energy_lines[1] - k_energy_lines[0])

	for i in range(1, len(k_energy_lines)):

		tmp = []
		for ii in range(k_energy_lines[i-1]+3, k_energy_lines[i]):
			lines = linecache.getline(logfile, ii).split()
			for iii in lines:
				tmp.append(float(iii))
		k_energy.append(tmp)

	if len(k_coordinates_c) == 2*len(k_energy):
		k_coordinates_c = k_coordinates_c[:91]

	if len(k_coordinates_c) == len(k_energy):
		print('The number of k-points = ', len(k_energy))
		print('The number of bands = ', len(k_energy[0]))
	else:
		print(len(k_energy))
		print(len(k_coordinates_c))
		print('k_coordinates) != len(k_energy)')
		print('Exitting...')
		exit()

	lattice_r = np.array(lattice_r)
	print('lattice_r (reciprocal axes)=', '\n', lattice_r)

	k_coordinates_c = np.array(k_coordinates_c); k_energy = np.array(k_energy)

	k_coordinates_d = np.dot(k_coordinates_c, np.linalg.inv(lattice_r))

	# write k_distance
	k_distance = [0]; sum_ = 0; k_energy = k_energy - fermi * au
	for i in range(1, len(k_energy)):
		k_pre = k_coordinates_c[i-1]
		k_aft = k_coordinates_c[i]
		dis = np.linalg.norm(k_aft-k_pre)
		sum_ += dis
		k_distance.append(round(sum_,6))

	# write BAND.dat 
	with open('BAND.dat', 'w+') as f:

		f.writelines('#K-Path(1/A) Energy-Level(eV)'+'\n')
		f.writelines('# NKPTS & NBANDS:  %s  %s' %(len(k_energy), len(k_energy[0]))+'\n')

		for i in range(len(k_energy[0])): # bands
			f.writelines('# Band-Index     %s' %(i+1)+'\n')
			if i%2 == 0: # positive sequence, kpoints
				for ii in range(len(k_energy)):
					f.writelines('    ' + '%0.4f' %k_distance[ii] + '    ' + '%0.4f' %k_energy[ii, i]+'\n')
				f.writelines('\n')

			if i%2 == 1: # negative sequence, kpoints
				for ii in range(len(k_energy)-1,-1,-1):
					f.writelines('    ' + '%0.4f' %k_distance[ii] + '    ' + '%0.4f' %k_energy[ii, i]+'\n')
				f.writelines('\n')

	# write KLABELS
	with open('KLABELS', 'w+') as f:
		f.writelines('K-Label    K-Coordinate in band-structure plots '+'\n')
		f.writelines(K_lable[0]+'                  '+'0'+'\n')
		for i in range(1, len(K_lable)):
			f.writelines(K_lable[i]+'                  ')
			for ii in range(1,len(k_energy)):
				dis = np.linalg.norm(k_coordinates_d[ii,:]-K_coord[i])
				if abs(dis) < 10**-3:
					f.writelines(str(k_distance[ii])+'\n')
					
		f.writelines('\n'+ 'Give the label for each high symmetry point in KPOINTS (KPATH.in) file. Otherwise, they will be identified as \'Undefined\' in KLABELS file')

if Getband_fun['TorF']:
	print('**********Running Getband!!!**********')
	if os.path.isfile(Getband_fun['log_file']):
		Getband(logfile = Getband_fun['log_file'], fermi = Getband_fun['fermi'],
			    K_lable = Getband_fun['K_lable'], K_coord = Getband_fun['K_coord'])
	else:
		print('There is no %s!!!' % Getband_fun['log_file'])

######################################
def Getwfn(inputs, atom_numb, atom_kind):
	
	# read positions from inputs file
	with open(inputs, 'r') as f:
		data = f.readlines()
		for lines in data:
			if 'CELL_PARAMETERS angstrom' in lines:
				cell_index = data.index(lines)
				print('cell_index =',cell_index)
			if 'nat' in lines:
				total_atom_num = int(lines.split()[2])
				print('total_atom_num =', total_atom_num)
			if 'ATOMIC_POSITIONS crystal' in lines:
				position_index = data.index(lines)
				print('position_index =', position_index)
	
	lattice = [linecache.getline(inputs, i+cell_index+2) for i in range(3)]
	
	position = [linecache.getline(inputs, i+position_index+2)[2:] for i in range(total_atom_num)]
	
	# write POSCAR
	with open('POSCAR.vasp', 'w+') as f:
		f.writelines('POSCAR'+'\n'+'1.0'+'\n')
		for i in lattice:
			f.writelines(i)
		
		f.writelines(atom_kind+'\n'+atom_numb+'\n'+'Direct'+'\n')
		
		for i in position:
			f.writelines(i)

if Getwfn_fun['TorF']:
	print('**********Running Getband!!!**********')
	if os.path.isfile(Getwfn_fun['inputs']):
		Getwfn(inputs = Getwfn_fun['inputs'],
			   atom_kind = Getwfn_fun['atom_kind'],
			   atom_numb = Getwfn_fun['atom_numb'])
	else:
		print('There is no %s or %s!' %(Getwfn_fun['inputs']))



