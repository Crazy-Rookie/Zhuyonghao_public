#! /usr/bin/python3
# -*- conding=UTF-8 -*-
'''
1. I can correctly calculate distance between two atoms.
2. Support only Direct type of POSCAR. 
3. author: yonghao_zhu@163.com
'''
##################################
import numpy as np
import linecache, os
##################################
'''
1. input files: POSCAR with Direct mode
2. Atoms: start 1
'''
Distancefun          = {}
Distancefun['TorF']  = 'T' #T
Distancefun['Atom1'] = 4
Distancefun['Atom2'] = 6
##################################
def ReadPOSCAR(file='POSCAR',print_='T'):
    '''
    input files: POSCAR
    '''
    atoms_kind = linecache.getline('POSCAR',6).split()
    atoms_num = [int(i) for i in linecache.getline('POSCAR',7).split()]
    #test#print(atoms_kind,atoms_num)
    if print_ == 'T':
        print('*******POSCAR*******')
        print('',atoms_kind,'\n',atoms_num)
        for i in range(len(atoms_kind)):
            if i == 0:
                print('',atoms_kind[i],'-->','1 -',atoms_num[i])
            if i > 0:
                print('',atoms_kind[i],'-->',sum(atoms_num[:i])+1,'-',sum(atoms_num[:i])+atoms_num[i])
        print('********************')
    comments = [linecache.getline(file,i).rstrip('\n') for i in range(9)]
    del comments[0]
    pos_np = np.loadtxt(file,comments=comments)
    return atoms_kind,atoms_num,pos_np
    #test#print(pos_np.shape)

def ModificationCoordinate(pos):
    '''
    Support Only Direct Type
    pos.shape = [NumAtoms,3]
    '''
    pos_new = np.zeros_like(pos)
    for i in range(pos.shape[0]):
        for m in range(3): 
            if pos[i,m] > 1:
                pos_new[i,m] = pos[i,m] - 1
            if pos[i,m] < 0:
                pos_new[i,m] = pos[i,m] + 1
            if pos[i,m] >= 0 and pos[i,m] <= 1:
                pos_new[i,m] = pos[i,m]

    return pos_new

def Distance(Atom1=1,Atom2=2):

    #read POSCAR and lattice
    pos_old = ReadPOSCAR(file='POSCAR',print_='F')[2]
    latt = np.array([[float(m) for m in linecache.getline('POSCAR',i).split()] for i in range(3,6)])
    #test#print(pos_old.shape)
    pos_new = ModificationCoordinate(pos=pos_old)
    #test#print(pos_new.shape)

    #make supercell (3,3,3)
    supercell_ = []
    for m in range(-1,2):
        for n in range(-1,2):
            for q in range(-1,2):
                a = []
                a.append(m)
                a.append(n)
                a.append(q)
                supercell_.append(a)
    supercell_ = np.array(supercell_).reshape(27,3)
  
    pos_Atom1 = (pos_new[Atom1-1] + supercell_).dot(latt)
    pos_Atom2 = pos_new[Atom2-1].dot(latt)

    distance_ = np.linalg.norm(pos_Atom1 - pos_Atom2,axis = 1)# for i in range(pos_Atom1.shape[0]) 
    distance = np.min(distance_)
    
    print('###################################')
    print('Atoms: start 1')
    print('Distance Atom%s-Atom%s : %.04f' %(Atom1,Atom2,distance),'Å')
    print('###################################')
    
if Distancefun['TorF'] == 'T' and os.path.isfile('POSCAR'):
    Distance(Atom1=Distancefun['Atom1'],Atom2=Distancefun['Atom2'])

