#! /usr/bin/python3
# -*- conding=UTF-8 -*-

'''
1. I can correctly calculate Standard Deviations(Å).
2. Support only Direct type of XDATCAR. 
3. author: yonghao_zhu@163.com
'''

#############################
import CorrectDistance as CD
import numpy as np
import math
#############################

#############################################
SD_fun = {}
SD_fun['TorF'] = True
#time_ini/time_fin: including, start 1
#time_ini is standard. 
SD_fun['time_ini']  = 1
SD_fun['time_fin']  = 10
SD_fun['atom_list'] = [37,    37 ,  148]
SD_fun['save_avg']  = False  # True or False
SD_fun['save_time'] = False  # True or False
#############################################

def PrintOut_avg(atom_list = [37, 148, 37], sd_atoms = np.zeros([1,1])):

	for i in range(len(atom_list)):
		if i == 0:
			print(atom_list[i],'--','1 -',atom_list[i])
			tmp = sd_atoms[: atom_list[i]]
			avg = sum(tmp) / atom_list[i]
			print(avg)
			print('=========================')

		if i > 0:
			print(atom_list[i],'---',sum(atom_list[:i])+1,'-',sum(atom_list[:i])+atom_list[i])
			tmp = sd_atoms[sum(atom_list[:i]): sum(atom_list[:i])+atom_list[i]]
			avg = sum(tmp) / atom_list[i]
			print(avg)
			print('=========================')

def SaveOut_time(atom_list = [37, 148, 37], sd_atoms_time = np.zeros([1,1])):

	print('sd_atoms_time.shape=[atoms, times]=',sd_atoms_time.shape)

	for i in range(len(atom_list)):
		name = '%s.txt' % (i+1)
		if i == 0:
			tmp = sd_atoms_time[: atom_list[i], :]
		if i > 0:
			tmp = sd_atoms_time[sum(atom_list[:i]): sum(atom_list[:i])+atom_list[i], :]

		avg = np.mean(tmp, axis=0)
		np.savetxt(name, avg)


def StandardDeviations(
	file = 'XDATCAR', 
	time_ini = 1, time_fin = 1000, 
	save_avg = False, save_time = False, 
	atom_list = [37, 148, 37]
	):
	
	'''
	time_ini/time_fin: including, start 1
	time_ini is standard time. 
	'''

	xdatcar, lattice, atoms = CD.ReadXDATCAR(time_ini = time_ini, time_fin = time_fin,)
	#xdatcar.shape = (times, atoms, 3)

	distance_1_2,pos_new = np.zeros([time_fin-time_ini+1, atoms, 3]), np.zeros([time_fin-time_ini+1, atoms, 3])
	#pos_new.shape = (times, atoms, 3) 

	sd_atoms = []; sd_atoms_time = np.zeros([atoms, time_fin-time_ini+1])

	for atom in range(atoms):

		pos_atom = np.array([xdatcar[time,atom,:] for time in range(time_fin - time_ini + 1)])
		#pos_atom.shape = (times, atom)

		pos_new[0, atom, :] = pos_atom[0, :]
		
		for time in range(1, time_fin - time_ini + 1):
		 	pos_tmp = np.stack([pos_atom[time], pos_atom[0]])
		 	#return Atom1(pos_atom[time])
		 	pos_new[time, atom, :] = CD.Distance(pos=pos_tmp, print_='F', latt=lattice,
		 		                          Atom1=1, Atom2=2)[1]

		pos_time_atom = pos_new[:, atom, :]
		#pos_time_atom.shape = (times, 3)
		
		avg_pos = np.mean(pos_time_atom, axis=0)
		#avg_pos.shape = (1, 3)

		distance_ = (pos_time_atom - avg_pos).dot(lattice)
		distance = np.linalg.norm(distance_, axis=1)
		distance_2 = distance**2
		#distance.shape = (times, )

		# print(distance.shape)
		# print('------------------')

		sd_atoms.append(math.sqrt(np.mean(distance_2)))

		if save_avg:
			with open('SD_atoms.txt', 'a+') as f:
				f.writelines(str(np.mean(distance))+'\n')

		if save_time:
			sd_atoms_time[atom, :] = distance

	if save_time and len(atom_list) != 0:
		SaveOut_time(atom_list = atom_list, sd_atoms_time = sd_atoms_time)

	if len(atom_list) != 0:
		PrintOut_avg(atom_list = atom_list, sd_atoms = sd_atoms)

def main():
	if SD_fun['TorF']:
		StandardDeviations(file = 'XDATCAR',
						   time_ini = SD_fun['time_ini'], time_fin = SD_fun['time_fin'],
						   save_avg = SD_fun['save_avg'], save_time = SD_fun['save_time'],
					       atom_list = SD_fun['atom_list'])


if __name__ == '__main__':
	main()

