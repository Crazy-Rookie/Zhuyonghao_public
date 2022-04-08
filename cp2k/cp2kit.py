#! /usr/bin/python3
# --* coding UTF-8 *--
'''
    0. pre-processing script for cp2k
    1. author: yonghao_zhu@163.com
'''
########################################################
import sys, os, math
import linecache
import numpy as np
import pandas as pd
########################################################
#                                                      #
#                geo_opt and cell_opt                  #
#                                                      #
#------------------------------------------------------#
##task 11: SUBSYS--COORD(VASPtoXYZ)b
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
VASPtoXYZ_fun['TorF'] = 'F'
#------------------------------------------------------#
##task 12: get xyz file from GeoOpt-pos-1.xyz
GetFinXTZ_fun = {}
'''
    input file:  GeoOpt-pos-1.xyz or CellOpt-pos-1.xyz
    output file: CONTCAR.xyz
'''
GetFinXTZ_fun['TorF'] = 'F'
########################################################
########################################################
#                                                      #
#                         band                         #
#                                                      #
#------------------------------------------------------#
##task 3: rewrite band.bs
'''
    input file:     band.bs
    output file:    band_plot.dat
    temperature:    0K
    unit of energy: 1au = 27.2113838565563 eV
'''
ReadBand_fun            = {}
ReadBand_fun['TorF']    = 'F'
#get fermi level in charge.out --> grep Fermi charge.out
ReadBand_fun['E_fermi'] = 7.084541 / 2.72113838565563E+01 / 2 #[Ry]
#cell_A, cell_B and cell_C: angstrom
ReadBand_fun['cell_A']  = [[9.2319517136, 0.0000000000, 0.0000000000]]
ReadBand_fun['cell_B']  = [[0.0000000000, 9.2319507599, 0.0000000000]]
ReadBand_fun['cell_C']  = [[0.0000000000, 0.0000000000, 8.8956279755]]
########################################################
########################################################
#                                                      #
#                        MD                            #
#                                                      #
#------------------------------------------------------#
##task 21-MD(NVT): Time-Temp(K) and Time-TotalEne(au)
'''
    input file:     GET_MD-1.ener
    output file:    Time-Temp.dat and Time-Ene.dat 
    unit of energy: 1au = 27.2113838565563eV
'''
AnalysisMD_21 = {}
AnalysisMD_21['TorF'] = 'F'
#------------------------------------------------------#
##task 22-MD(NVT): Time-BandEne(eV)
'''
    input file:   md.out
    output file:  Time-BandEnergy.dat 
    unit of energy: 1au = 27.2113838565563eV
'''
AnalysisMD_22          = {}
AnalysisMD_22['TorF']  = 'F'
AnalysisMD_22['POTIM'] = 1 #[fs]
AnalysisMD_22['NHOMO'] = 2 #include HOMO
AnalysisMD_22['NLUMO'] = 2 #include LUMO
#------------------------------------------------------#
##task 23-MD(NVT): split trajectory
AnalysisMD_23          = {}
'''
    input file:   GET_MD-pos-1.xyz
    output file:  ./pfiles/p****.xyz    
'''
AnalysisMD_23['TorF']   = 'F'
AnalysisMD_23['length'] = 4
#------------------------------------------------------#
##task 24-MD(NVT): convert cp2k trajectory into vasp trajectory
##                 GET_MD-pos-1.xyz to XDATCAR
AnalysisMD_24           = {}
AnalysisMD_24['TorF']   = 'F'
'''
    input file:   ./pfiles/p****.xyz
    output file:  ./vasp_pfiles/p****
'''
AnalysisMD_24['cell']   = ['45.9291000366 0.0000000000 0.0000000000',
                           '0.0000000000 32.4765014648 0.0000000000',
                           '0.0000000000 0.0000000000 32.1833000183']
#------------------------------------------------------#
##task 25: convert XDATCAR to MD_traj.xyz
AnalysisMD_25           = {}
AnalysisMD_25['TorF']   = 'F'
########################################################
########################################################
#                                                      #
#                   FUNCTIONS RUN                      #
#                                                      #
########################################################
########################################################
#                                                      #
#                geo_opt and cell_opt                  #
#                                                      #
#------------------------------------------------------#
#task 11: SUBSYS--COORD(VASPtoXYZ)
def VASPtoXYZ(file_VASP='POSCAR.vasp',
    file_XYZ='POSCAR.xyz'):
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
    pos_np = np.loadtxt(file_VASP,comments=comments)
    print('POSCAR shape:',pos_np.shape)
    
    if 'D' in linecache.getline(file_VASP,8): #Direct
        print('Direct......')
        latt_tmp = [linecache.getline(file_VASP,i).rstrip('\n').split()
            for i in range(3,6)]
        latt = np.array([[float(i[0]),float(i[1]),float(i[2])] for i in latt_tmp])
        xyz = np.dot(pos_np,latt)
        print(xyz)

    if 'C' in linecache.getline(file_VASP,8): #Cartesian
        print('Cartesian......')
        xyz = pos_np

    #write xyz file
    print('********************')
    for i in range(len(atoms_kind)):
        if i == 0:
            print('',atoms_kind[i],'-->','1 -',atoms_num[i])
            for ii in range(0,atoms_num[i]):
                xyz_file.writelines(atoms_kind[i]+'    '+
                    '%.9f' % xyz[ii][0]+'  '+
                    '%.9f' % xyz[ii][1]+'  '+
                    '%.9f' % xyz[ii][2]+'  '+'\n')
        if i > 0:
            print('',atoms_kind[i],'-->',sum(atoms_num[:i])+1,'-',sum(atoms_num[:i])+atoms_num[i])
            for ii in range(sum(atoms_num[:i]),sum(atoms_num[:i])+atoms_num[i]):
                xyz_file.writelines(atoms_kind[i]+'    '+
                    '%.9f' % xyz[ii][0]+'  '+
                    '%.9f' % xyz[ii][1]+'  '+
                    '%.9f' % xyz[ii][2]+'\n')
    print('********************')

if VASPtoXYZ_fun['TorF'] == 'T':
    if os.path.isfile('POSCAR.vasp'):
        VASPtoXYZ(file_VASP='POSCAR.vasp', file_XYZ='POSCAR.xyz')
    else:
        print('No POSCAR.vasp!!!')
#------------------------------------------------------#
#task 12: get xyz file from GeoOpt-pos-1.xyz
def GetFinXTZ(file_inp='GeoOpt-pos-1.xyz',
    file_out='CONTCAR.xyz'):
    file_i = open(file_inp,'r')
    atoms_num = int(linecache.getline(file_inp,1))
    print('atoms_num:',atoms_num)
    pos = []
    for line in file_i:
        lines = line.split()
        if len(lines) == 4:
            pos.append(line)
    file_i.close()
    pos = np.array(pos).reshape(-1,atoms_num)
    print('pos.shape(steps,atoms):',pos.shape)
    file_o = open(file_out,'w+')
    file_o.writelines(str(atoms_num)+'\n'+'\n')
    for i in pos[-1,:]:
        file_o.writelines(i)
    file_o.close()

if GetFinXTZ_fun['TorF'] == 'T':
    if os.path.isfile('GeoOpt-pos-1.xyz'):
        GetFinXTZ(file_inp='GeoOpt-pos-1.xyz',
            file_out='CONTCAR.xyz')
    elif os.path.isfile('CellOpt-pos-1.xyz'):
        GetFinXTZ(file_inp='CellOpt-pos-1.xyz',
            file_out='CONTCAR.xyz')  
    else:
        print('No GeoOpt-pos-1.xyz!!!')
########################################################
########################################################
#                                                      #
#                        MD                            #
#                                                      #
#------------------------------------------------------#
##task 2-MD(NVT): Time-Temp(K), Time-TotalEne(eV)
##Time-BandEnergy(eV)
def AnalysisMD(task=21, POTIM = 1, length = 1,
    file_inp='test.dat', file_temp='test.dat', file_ene='test.dat',
    NHOMO = 1, NLUMO = 1, 
    cell = ['1 0 0', '0 1 0', '0 0 1']):

    au = 2.72113838565563E+01
    
    #task 21
    if task == 21:
        np_inp = np.loadtxt(file_inp)
        
        times = np_inp[:,1]

        temp = np.zeros([times.shape[0],2])

        temp[:,0] = np_inp[:,1]; temp[:,1] = np_inp[:,3]
        np.savetxt(file_temp,temp,fmt='%.04f')

        ene = np.zeros([times.shape[0],2])
        ene[:,0] = np_inp[:,1]; ene[:,1] = np_inp[:,4]
        np.savetxt(file_ene,ene,fmt='%.08f')
    
    #task 22
    if task == 22:

        lines_lumo, lines_homo = [],[]
        if NLUMO%4 == 0:
            lines_lumo.append(4)
            lines_lumo.append(int(NLUMO/4))
        else:
            lines_lumo.append(NLUMO%4)
            lines_lumo.append(int(NLUMO/4) + 1)

        if NHOMO%4 == 0:
            lines_homo.append(4)
            lines_homo.append(int(NHOMO/4))
        else:
            lines_homo.append(NHOMO%4)
            lines_homo.append(int(NHOMO/4) + 1)

        md_out = open(file_inp, 'r', encoding='utf-8')
        fermi,bandene = [],[] #eV
        for num,line in enumerate(md_out) :
            if 'Fermi Energy [eV]' in line:
                #fermi energy [eV]
                fermi.append(float(line.split()[-1]))

                num = num + 1; tmp_be = []

                #get homos
                try:
                    
                    if lines_homo[1] == 1:
                        homo_ = linecache.getline(file_inp,num-1).split()
                        for ene in homo_[-1*lines_homo[0]:]:
                            
                            tmp_be.append(float(ene) * au) #au -> eV
                    else:
                        homo_ = linecache.getline(file_inp,num-lines_homo[1]).split()
                        for ene in homo_[-1*lines_homo[0]:]:
                            tmp_be.append(float(ene) * au) #au -> eV

                        for i in range(lines_homo[1]-1):
                            ii = lines_homo[1] - 1 - i
                            homo_ = linecache.getline(file_inp,num-ii).split()
                            for ene in homo_:
                                tmp_be.append(float(ene) * au) #au -> eV                            
                except:
                    print('ERROR, CHECK md.out!')

                #get lumos
                try:
                    if lines_lumo[1] == 1:
                        lumo_ = linecache.getline(file_inp,num+5).split()
                        for ene in lumo_[:lines_lumo[0]]:
                            
                            tmp_be.append(float(ene) * au) #au -> eV
                    else:

                        for i in range(lines_lumo[1]-1):

                            lumo_ = linecache.getline(file_inp,num+5+i).split()
                            for ene in lumo_:
                                tmp_be.append(float(ene) * au) #au -> eV  

                        lumo_ = linecache.getline(file_inp,num+5+lines_lumo[1]-1).split()
                        for ene in lumo_[:lines_lumo[0]]:
                            tmp_be.append(float(ene) * au) #au -> eV
                except:
                    print('ERROR, CHECK md.out!')
                
                if len(tmp_be) == NLUMO+NHOMO:
                    bandene.append(tmp_be)
                print(len(tmp_be) == NLUMO+NHOMO)
        
        md_out.close()

        band_name_h = ['HOMO-%s' % i for i in range(NHOMO)]
        band_name_l = ['LUMO+%s' % i for i in range(NLUMO)] 
        print(len(bandene))
        time = np.array([i for i in range(len(bandene))]) * POTIM
        bandene = np.array(bandene)
        bandene = np.insert(bandene, 0, values=time, axis=1)
        np.savetxt(file_ene, bandene, fmt='%0.5f')   

        with open(file_ene,'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('#Time(fs)')
            
            for i in range(1,len(band_name_h)+1):
                f.write('   '+band_name_h[-1*i])

            for i in band_name_l:
                f.write('   '+i)

            f.write('\n'+content)

    #task 23
    if task == 23:
        file_i = open(file_inp,'r')
        atoms_num = int(linecache.getline(file_inp,1))
        print('atoms_num:',atoms_num)
        pos = []
        for line in file_i:
            lines = line.split()
            if len(lines) == 4:
                pos.append(line)
        file_i.close()
        pos = np.array(pos).reshape(-1,atoms_num)
        print('pos.shape(steps,atoms):',pos.shape)

        for i in range(pos.shape[0]):
            try:
                xx = './pfiles/p%03d.xyz'.replace('3', str(length))
                file_name = xx % i
                with open(file_name,'w+') as f:
                    for ii in pos[i]:
                        f.write(ii)
            except:
                print('ERROR, KEEP ./pfiles/!!!')

    #task 24
    if task == 24:
        pd.set_option('display.max_columns',None)

        pfiles_name = [i for i in os.listdir('./pfiles/') if 'p' in i]

        for pfiles in pfiles_name:
            p_df = pd.read_csv('./pfiles/' + pfiles,sep="\s+",header=None,names=['a','b','c','d'])

            a_value = p_df['a'].value_counts()
            a_counts = pd.DataFrame(a_value).reset_index().rename(columns={'index':'name','a':'counts'})

            pos_ = np.zeros([len(p_df),3])

            atom_types = a_counts.values

            index_ = 0
            for i in atom_types[:,0]:
                atom_i = p_df[p_df['a']==i][['b','c','d']].values
                atom_num = atom_i.shape[0]
                pos_[index_: index_ + atom_num, :] = atom_i
                index_ += atom_num

            if not os.path.isfile('./vasp_pfiles/' + pfiles.replace('.xyz','')):
                np.savetxt('./vasp_pfiles/' + pfiles.replace('.xyz',''), pos_, fmt = '%08f')

            project_name = 'vasp pfiles'

            with open('./vasp_pfiles/' + pfiles.replace('.xyz',''), 'r+') as f:
                content = f.read()        
                f.seek(0, 0)
                f.write(project_name + '\n' + '1' + '\n')
                for i in cell:
                    f.write(i + '\n')
                for i in range(atom_types.shape[0]):
                    f.write(atom_types[i,0] + '  ')
                f.write('\n')
                for i in range(atom_types.shape[0]):
                    f.write(str(atom_types[i, 1]) + '  ')
                f.write('\n' + 'C' + '\n' + content)

    #task 25
    if task == 25:
        seventh_line = linecache.getline(file_inp,7).split()
        atoms_kind = [int(i) for i in seventh_line]
        atoms = sum(atoms_kind)
        comments = ['Direct configuration=']
        for i in range(7):
            comments.append(linecache.getline(file_inp,i+1))
        xdatcar = np.loadtxt(file_inp,comments=comments).reshape(-1,atoms,3)
        print('XDAtCAR.shape =', xdatcar.shape)

        lattice = np.zeros([3,3])
        for i in range(3):
            tmp = [float(j) for j in comments[3+i].split()]
            lattice[i,:] = tmp

        xdatcar = np.dot(xdatcar, lattice)

        elements = linecache.getline(file_inp,6).split()
        
        ele_list = []
        for i in range(len(atoms_kind)):
            num = atoms_kind[i]
            tmp = [elements[i]] * num
            for j in tmp:
                ele_list.append(j)
        
        with open('MD_traj.xyz', 'w+') as f:

            tmp_0 = ' i =        x, time =        x.000, E =       -100'

            for time in range(xdatcar.shape[0]):
                f.writelines('  ' + str(atoms) + '\n')
                tmp_1 = tmp_0.replace('x', str(time))
                f.writelines(tmp_1 + '\n')
                pos = xdatcar[time]
                for i in range(atoms):
                    tmp_2 = ' ' + ele_list[i] + '         '
                    tmp_2 = tmp_2 + '%.10f' %pos[i, 0] + '         ' + '%.10f' %pos[i, 1] + '         ' + '%.10f' %pos[i, 2]
                    f.writelines(tmp_2 + '\n')
##run
#task-21
if AnalysisMD_21['TorF'] == 'T':
    if os.path.isfile('GET_MD-1.ener'):
        AnalysisMD(task = 21,
            file_inp = 'GET_MD-1.ener',
            file_temp = 'Time-Temp.dat',
            file_ene = 'Time-TotalEne.dat')
    else:
        print('No GET_MD-1.ener!!!')
#task-22
if AnalysisMD_22['TorF'] == 'T':
    if os.path.isfile('md.out'):
        AnalysisMD(task = 22,
            file_inp = 'md.out', POTIM = AnalysisMD_22['POTIM'],
            file_ene = 'Time-BandEnergy.dat',
            NHOMO = AnalysisMD_22['NHOMO'],
            NLUMO = AnalysisMD_22['NLUMO'])
    else:
        print('No files!!!')
#task-23
if AnalysisMD_23['TorF'] == 'T':
    if os.path.isfile('GET_MD-pos-1.xyz'):
        AnalysisMD(task=23, length=AnalysisMD_23['length'],
        file_inp='GET_MD-pos-1.xyz')
    else:
        print('No GET_MD-pos-1.xyz!!!')
#task-24
if AnalysisMD_24['TorF'] == 'T':
    AnalysisMD(task=24, cell = AnalysisMD_24['cell'])
#task-25
if AnalysisMD_25['TorF'] == 'T':
    AnalysisMD(task=25, file_inp='XDATCAR')
########################################################
########################################################
#                                                      #
#                         band                         #
#                                                      #
#------------------------------------------------------#
##task 3: rewrite band.bs
def ReadBand(file_in='band.bs',file_out='BAND.dat',fermi=1,
            cell_A=[1,0,0],
            cell_B=[0,1,0],
            cell_C=[0,0,1]):

    cell_A = np.array(cell_A)
    cell_B = np.array(cell_B)
    cell_C = np.array(cell_C)
    cell_v = np.dot(cell_A,np.cross(cell_B,cell_C).T)
    print('volume of cell= ',round(cell_v[0][0],6),'angstrom^3')

    cell_A_reci = np.abs(np.cross(cell_B,cell_C)/cell_v)
    cell_B_reci = np.abs(np.cross(cell_A,cell_C)/cell_v)
    cell_C_reci = np.abs(np.cross(cell_B,cell_A)/cell_v)
    print('-------------------------------------')
    print('reciprocal lattice vectors(like vasp):')
    cell_reci = np.vstack((cell_A_reci,cell_B_reci,cell_C_reci))
    print(cell_reci)
    print('-------------------------------------')
    
    comments = []
    for i in range(1,100):
        line = linecache.getline(file_in,i)
        if 'Nr.' not in line:
            comments.append(line)
        else:
            break

    comments.append(linecache.getline(file_in,2+len(comments)))
    comments.append('Nr.')

    nkpoints = int(comments[0].split()[-1])
    nbands = int(comments[-2])
    print('nbands=',nbands)
    print('nkpoints=',nkpoints)
    
    au = 27.2113838565563
    if nbands % 4 == 0:
        band_pri = np.loadtxt(file_in,comments=comments).reshape(nkpoints,nbands) - fermi*au/2
    else:
        with open(file_in) as f:
            band_pri = []
            lines = f.readlines()
            for ii in range(len(lines)):
                if 'Nr.' in lines[ii]:
                    band_tmp = []
                    for band_line in range(int(nbands/4) + 1):
                        for iii in lines[ii+2+band_line].split():
                            band_tmp.append(float(iii))
                    band_pri.append(band_tmp)
            band_pri = np.array([band_pri]).reshape(nkpoints,nbands) - fermi*au/2

    kpoint_1 = []
    with open(file_in) as f:
        for line in f:
            if 'Nr.' in line:
                tmp = [float(line.split()[5]),
                    float(line.split()[6]),
                    float(line.split()[7])]
                kpoint_1.append(tmp)
    kpoint_1 = np.dot(np.array(kpoint_1).reshape(nkpoints,3),cell_reci)
    
    print('kpoint_1.shape=',kpoint_1.shape)

    kpoint_2 = [0]
    for i in range(0,kpoint_1.shape[0]-1):
        kpoint_1_1 = kpoint_1[i]
        kpoint_1_2 = kpoint_1[i+1]
        tmp_0 = np.linalg.norm(kpoint_1_1 - kpoint_1_2)*2*math.pi #like vasp
        tmp_1 = kpoint_2[i] + tmp_0
        kpoint_2.append(round(tmp_1,6))

    file_o = open(file_out,'w+')
    file_o.writelines('#K-Path      Energy-Level'+'\n'+'# NKPTS & NBANDS:  '+
        str(nkpoints)+'  '+str(nbands)+'\n')

    print(len(kpoint_2))
    print(band_pri.shape)

    for iband in range(nbands):
        file_o.writelines('# Band-Index    '+str(iband)+'\n')
        for ikpt in range(nkpoints):
            file_o.writelines(str(kpoint_2[ikpt])+'    '+str(band_pri[ikpt,iband])+'\n')
        file_o.writelines('\n')
        
    file_o.close()

if ReadBand_fun['TorF'] == 'T':
    if os.path.isfile('band.bs'):
        ReadBand(file_in='band.bs',file_out='band_plot.dat',
            fermi=ReadBand_fun['E_fermi'],
            cell_A=ReadBand_fun['cell_A'],
            cell_B=ReadBand_fun['cell_B'],
            cell_C=ReadBand_fun['cell_C'])
    else:
        print('No band.bs file!!!')
########################################################
########################################################
########################################################


########################################################


########################################################
########################################################
#                                                      #
#                         unit                         #
#                                                      #
#------------------------------------------------------#
'''
    import sys
    sys.path.append('D:\\360\\anaconda3\\Lib\\site-packages')
    
    [u] -> [a.u.]                   1.82288848426455E+03
    [Angstrom] -> [Bohr] = [a.u.]   1.88972613288564E+00
    [a.u.] = [Bohr] -> [Angstrom]   5.29177208590000E-01
    [a.u.] -> [s]                   2.41888432650478E-17
    [a.u.] -> [fs]                  2.41888432650478E-02
    [a.u.] -> [J]                   4.35974393937059E-18
    [a.u.] -> [N]                   8.23872205491840E-08
    [a.u.] -> [K]                   3.15774647902944E+05
    [a.u.] -> [kJ/mol]              2.62549961709828E+03
    [a.u.] -> [kcal/mol]            6.27509468713739E+02
    [a.u.] -> [Pa]                  2.94210107994716E+13
    [a.u.] -> [bar]                 2.94210107994716E+08
    [a.u.] -> [atm]                 2.90362800883016E+08
    [a.u.] -> [eV]                  2.72113838565563E+01
    [a.u.] -> [Hz]                  6.57968392072181E+15
    [a.u.] -> [1/cm] (wave numbers) 2.19474631370540E+05
    [a.u./Bohr**2] -> [1/cm]        5.14048714338585E+03
'''