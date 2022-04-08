# /usr/bin/python3

'''
    The python3 script was written to process the MD and NAC data.
    Author email: yonghao_zhu@163.com
'''

###################################
import os
import re
import linecache
import math
import numpy as np
import sys
import getopt
###################################

###################################
#-task 1
def TotalEne(file='OSZICAR'):
    if os.path.isfile(file):
        oszicar = open(file,'r');time = 1
        movie = {}
        #read Energy
        for line in oszicar:
            line = line.split()
            if 'T=' in line:
                movie[time] = [float(line[8]),float(line[2])]
                time += 1
        ##test##print(movie[12])
        oszicar.close
        energy_file = open('TimeTemEn.dat','w+')
        energy_file.writelines('Time'+'        '+'Energy'+'        '+'Temperature'+'\n')
        time = [];energy = []
        for key in movie:
            time.append(key)
            energy.append(movie[key])
        ti = 0
        while ti < len(time):
            energy_file.writelines(str(time[ti])+'        '+str(energy[ti][0])+'        '+str(energy[ti][1])+'\n')
            ti += 1
        energy_file.close
    else:
        print('There is no %s!' % file)
###################################

###################################
#-task 21/22 -nkpt 1 -band_ini 176 -band_fin 177
def BandEne(task=21,nkpt=1,band_ini=176,band_fin=177,ispin=1):
    if task == 21:
        if os.path.isfile('OUTCAR'):
            outcar=open('OUTCAR')

            if ispin == 1: 
                band_file=open('band_energy.dat','w+')
                for n,line in enumerate(outcar):
                    if nkpt == 1: 
                        if " E-fermi" in line:
                            a1=linecache.getlines('OUTCAR')[n]
                            a2=a1.strip().split()[2]
                            a3=float(a2) #E-fermi
                            for i in range(band_ini,band_fin+1):
                                b2=linecache.getlines('OUTCAR')[n+4+i]
                                b3=b2.strip().split()[1]
                                b4=float(b3)
                                b5=str("%.4f"%(b4))
                                band_file.write(b5+' ')
                            band_file.write('\n')

                    if nkpt != 1:
                        if "NBANDS" in line:
                            nbands = int(line.split()[14])
                        if " E-fermi" in line:
                            a1=linecache.getlines('OUTCAR')[n]
                            a2=float(a1.strip().split()[2]) #E-fermi
                            for ikpt in range(nkpt):
                                for i in range(band_ini,band_fin+1):
                                    b2=linecache.getlines('OUTCAR')[(nbands+3)*ikpt+n+4+i]
                                    b3=b2.strip().split()[1]
                                    b4=float(b3)
                                    b5=str("%.4f"%(b4))
                                    band_file.write(b5+' ')
                            band_file.write('\n')
                            
                band_file.close()

            if ispin == 2:
                band_file_up = open('band_energy_UP.dat','w+')
                band_file_dw = open('band_energy_DW.dat','w+')
                for n,line in enumerate(outcar):
                    if "NBANDS" in line:
                        nbands = int(line.split()[14])
                    if nkpt == 1: 
                        if " E-fermi" in line:
                            a1=linecache.getlines('OUTCAR')[n]
                            a2=a1.strip().split()[2]
                            a3=float(a2) #E-fermi
                            for i in range(band_ini,band_fin+1):
                                #spin up
                                b2_up = linecache.getlines('OUTCAR')[n+6+i]
                                b3_up = b2_up.strip().split()[1]
                                b4_up = float(b3_up)
                                b5_up = str("%.4f"%(b4_up))
                                band_file_up.write(b5_up+' ')
                                #spin dw
                                b2_dw = linecache.getlines('OUTCAR')[n+11+i+nbands]
                                b3_dw = b2_dw.strip().split()[1]
                                b4_dw = float(b3_dw)
                                b5_dw = str("%.4f"%(b4_dw))
                                band_file_dw.write(b5_dw+' ')
                            band_file_up.write('\n')
                            band_file_dw.write('\n')

                    if nkpt != 1:
                        if " E-fermi" in line:
                            a1=linecache.getlines('OUTCAR')[n]
                            a2=float(a1.strip().split()[2]) #E-fermi
                            for ikpt in range(nkpt):
                                for i in range(band_ini,band_fin+1):
                                    #spin_up
                                    b2_up = linecache.getlines('OUTCAR')[(nbands+3)*ikpt+n+6+i]
                                    b3_up = b2_up.strip().split()[1]
                                    b4_up = float(b3_up)
                                    b5_up = str("%.4f"%(b4_up))
                                    band_file_up.write(b5_up+' ')
                                    #spin_dw
                                    b2_dw = linecache.getlines('OUTCAR')[(nbands+3)*ikpt+n+11+i]
                                    b3_dw = b2_dw.strip().split()[1]
                                    b4_dw = float(b3_dw)
                                    b5_dw = str("%.4f"%(b4_dw))
                                    band_file_up.write(b5_dw+' ')
                            band_file_up.write('\n')
                            band_file_dw.write('\n')

                band_file_up.close();band_file_dw.close()

            outcar.close()
        else:
            print('There is no OUTCAR!')
    if task == 22:
        if os.path.isfile('band_energy.dat'):
            ene = np.loadtxt('band_energy.dat')
            return ene.mean(axis=0)
        else:
            print('There is no %s' % 'band_energy.dat')
###################################

###################################
#-task 3 -time_ini 1000 -time_fin 2000 -save_ 0 -direct 1
def StandardDeviations(filename='XDATCAR',save_=0,time_ini=1000,time_fin=2000,direct=1):
    if os.path.isfile(filename):
        lattice = np.mat([[float(linecache.getline(filename,i+1).split()[0]),
            float(linecache.getline(filename,i+1).split()[1]),
            float(linecache.getline(filename,i+1).split()[2])] for i in range(2,5)])
        seventh_line = linecache.getline(filename,7).split()
        atoms = int(sum([float(i) for i in seventh_line]))
        comments = ['Direct configuration=']
        for i in range(7):
            comments.append(linecache.getline(filename,i+1))
        xdatcar_D = np.loadtxt(filename,comments=comments).reshape(-1,atoms,3)[time_ini:time_fin+1,:,:]
        xdatcar_C =np.array([np.dot(xdatcar_D[i],lattice) for i in range(time_fin-time_ini+1)])
        #STD = list(np.std(np.sqrt(np.sum(np.square(xdatcar_C-xdatcar_C.mean(axis=0)),axis=2)),axis=0))
        diff_dis = np.sum(np.square(xdatcar_C-xdatcar_C.mean(axis=0)),axis=2)
        STD = np.sqrt(np.sum(diff_dis,axis=0)/(time_fin-time_ini+1))
        temp = {};sum_atoms = 0
        if save_ == 1:
            print('The shape of STD: ',STD.shape,'\n','STD of Each atom:',STD)
            np.savetxt('STD_atoms.dat',STD)
            print('save STD_atoms.dat file.')
        for i in range(len(seventh_line)):
            key = 'elements%s' % i
            temp[key] = STD[sum_atoms:sum_atoms+int(seventh_line[i])]
            sum_atoms += int(seventh_line[i])
        out = []
        for i in temp:
            out.append(sum(temp[i])/len(temp[i]))
        print(out)
    else:
        print('There is no %s' %filename)
###################################

###################################
#-task 4 -md_times 6000
def split_xda(filename='XDATCAR'):
    if os.path.isfile(filename):
        seventh_line = linecache.getline(filename,7).split()
        atoms = int(sum([float(i) for i in seventh_line]))
        comments = ['Direct configuration='];lattice = []
        for i in range(7):
            comments.append(linecache.getline(filename,i+1))
            lattice.append(linecache.getline(filename,i+1))
        lattice.append('Direct \n')
        xdatcar = np.loadtxt(filename,comments=comments).reshape(-1,atoms,3)
        MD_times = xdatcar.shape[0]
        for i in range(MD_times):
            name = 'p' + '{0:0>4}'.format(i+1)
            print(i+1)
            file = open(name,'w+')
            file.writelines(j for j in lattice)
            file.writelines('     '+'{0:0<10}'.format(ii[0])+'    '+'{0:0<10}'.format(ii[1])+'    '+'{0:0<10}'.format(ii[2])+'\n' for ii in xdatcar[i])
    else:
        print('There is no %s' %filename)
###################################

###################################
#-task 51/52 -a1 75 -a2 76 -xdatcar 1 -car 0 time_ini 1 time_fin 2
def DistanceTwoAtoms(task=51,A1=75,A2=76,Car=1,time_ini=1,time_fin=2):
    file_dis = 'distance_%s_%s.dat' % (str(A1),str(A2))
    if task == 51:
        '''
            The distance between A1 and A2 changes with time.
            The POSCAR named 'p000X' obtained from pfiles.py
            NOTE: No '0.0000000 0.0000000 0.00000000' lines.
            Choice Car = 1 or 0.
        '''
        file = open(file_dis,'w+')
        for i in range(time_ini,time_fin+1):
            filename = 'p'+'{0:0>4}'.format(i)
            comments = [linecache.getline(filename,i+1)  for i in range(7)]        
            if Car == 1:
                comments.append('C')
                pos = np.loadtxt(filename,comments=comments)
                posA1 = pos[A1-1];posA2 = pos[A2-1]
                distance = np.sqrt(np.sum(np.square(posA2-posA1)))
                file.writelines(str(i+1)+'    '+str(distance)+'\n')
                print(filename)
            else:
                comments.append('D')
                lattice = np.mat([[float(linecache.getline(filename,i+1).split()[0]),
                    float(linecache.getline(filename,i+1).split()[1]),
                    float(linecache.getline(filename,i+1).split()[2])] for i in range(2,5)])
                pos = np.loadtxt(filename,comments=comments)
                posA1 = np.dot(pos[A1-1],lattice);posA2 = np.dot(pos[A2-1],lattice)
                distance = np.sqrt(np.sum(np.square(posA2-posA1)))
                file.writelines(str(i+1)+'    '+str(distance)+'\n')
                print(filename)            
        file.close()
    if task == 52:
        '''
            The distance between A1 and A2 changes with time.
            read the position in XDATCAR file.
            Car = 0
        '''
        if os.path.isfile('XDATCAR'):
            if Car == 0:
                seventh_line = linecache.getline('XDATCAR',7).split()
                lattice = np.mat([[float(linecache.getline('XDATCAR',i+1).split()[0]),
                    float(linecache.getline('XDATCAR',i+1).split()[1]),
                    float(linecache.getline('XDATCAR',i+1).split()[2])] for i in range(2,5)])
                comments = ['Direct configuration='];atoms = int(sum([float(i) for i in seventh_line]))
                for i in range(7):
                    comments.append(linecache.getline('XDATCAR',i+1))
                xdatcar = np.loadtxt('XDATCAR',comments=comments).reshape(-1,atoms,3)[time_ini:time_fin+1,:,:]
                posA1 = np.dot(xdatcar[:,A1-1,:],lattice);posA2 = np.dot(xdatcar[:,A2-1],lattice)
                distance = np.sqrt(np.sum(np.square(posA2-posA1),axis=1))
                np.savetxt(file_dis,distance,fmt='%.4e')
            else:
                print('Car should be 0 (Direct)!!!')
        else:
            print('There is no XDATCAR!')
    print('save %s file' % file_dis)
###################################

###################################
#-task 61/62/63/64/65 -time_fin 500 -time_fin 1000 -band_ini 52 -band_fin 53
def reduce_real(file='real',task_=60,band_ini=52,band_fin=53,time_ini=500,time_fin=1000):
    '''
        The function reads pristine real000x file after NAC calculate.
        And mkdir RedReal.
    '''
    if task_ == 61:#Reduce the pristine real files (NAC: real directory:61)
        os.popen('mkdir newmatrix')
        for i in range(time_ini,time_fin+1):
            filename = file+'{0:0>4}'.format(i)
            real_p = np.loadtxt(filename)
            real_r = real_p[band_ini:band_fin+1,2*band_ini:2*band_fin+2]
            print(filename)
            np.savetxt('./newmatrix/%s' % filename,real_r,fmt='%.4e')
        print('mkdir newmatrix, and reduce real files!')
    #----------------------------------------------------------
    if task_ == 62:#Rename the smaller real files (NAC: newmatrix directory:62)
        os.popen('mkdir -p %s_%s/res/' %(time_ini,time_fin))
        for i in range(time_ini,time_fin+1):
            filename_p = file+'{0:0>4}'.format(i)
            real_p = np.loadtxt(filename_p)
            filename_r = file+'{0:0>4}'.format(i-time_ini+1)
            print(filename_r)
            np.savetxt('%s_%s/res/%s' % (time_ini,time_fin,filename_r),real_p,fmt='%.4e')
        print('mkdir %s_%s/res' %(time_ini,time_fin))
    #----------------------------------------------------------
    if task_ == 63:#Reduce the band energy (NAC: band_energy.dat:63)
        ene_p = np.loadtxt('band_energy.dat')
        ene_r = ene_p[time_ini:time_fin+1]
        np.savetxt('reduce_energy',ene_r,fmt='%.4e')
    #----------------------------------------------------------
    if task_ == 64:#np.abs()absoluate value#Average NAC (NAC: smaller newmatrix/real:64)
        filename_1 = 'real'+'{0:0>4}'.format(time_ini)
        real_sum = np.abs(np.loadtxt(filename_1))
        for i in range(time_ini+1,time_fin+1):
            filename = file+'{0:0>4}'.format(i)
            real_sum += np.abs(np.loadtxt(filename))
        np.savetxt('average_NAC.dat',real_sum/(time_fin-time_ini+1),fmt='%.4e')
        print(real_sum/(time_fin-time_ini+1))
        print('save average.dat file!')
    #----------------------------------------------------------
    if task_ == 65:#time_num=time_max#np.abs()absoluate value#Average NAC (NAC: pristine real:65)
        filename_1 = 'real'+'{0:0>4}'.format(time_ini)
        real_sum = np.abs(np.loadtxt(filename_1)[band_ini:band_fin+1,2*band_ini:2*band_fin+2])
        for i in range(time_ini+1,time_fin+1):
            filename = file+'{0:0>4}'.format(i)
            real_sum += np.abs(np.loadtxt(filename)[band_ini:band_fin+1,2*band_ini:2*band_fin+2])
        np.savetxt('average_NAC.dat',real_sum/(time_fin-time_ini+1),fmt='%.4e')
        print(real_sum/(time_fin-time_ini+1))
        print('save average.dat file!')
    #----------------------------------------------------------
    if task_ == 66:#NAC-Time
        nac_time = open('NAC_Time.dat','w+')
        for i in range(time_ini,time_fin+1):
            filename = file+'{0:0>4}'.format(i)
            nac = str(np.loadtxt(filename)[band_ini,2*band_fin])
            print(band_ini,band_fin)
            print(str(np.loadtxt(filename)[band_ini,2*band_fin]))
            nac_time.writelines(str(i)+'    '+nac+'\n')
        nac_time.close()
###################################

###################################
#-task 7
def VaspComPyxaidEneNac(real_path='./res/',energy_file='./energy',
    name_num=4):
    '''
    1.combine the energy and nac files to Ham_*_im and Ham_*_re files
        for namd running in pyxaid, referring to combine.py file. 
    2.input files: ./real/real**** (full size) and ./energy (energy-time)
    3.progrom: re im re im --> re re + im im --> 0          re/(-1*au) (Ham_im_) + energy/au   im/(-1*au) (Ham_re_)
               re im re im     re re   im im     re/(-1*au) 0                      im/(-1*au)) energy/au
    4.1au = 13.60569253eV
    '''
    au = 13.60569253
    #get the number of real files and energy 
    time_real = len([i for i in os.listdir(real_path) if 'real' in i])
    energy = np.loadtxt(energy_file)/(au)
    time_energy = energy.shape[0]
    num_orbitals = energy.shape[1]
    #test#print('time_real=',time_real,'\ntime_energy=',time_energy,'\nnum_orbitals=',num_orbitals)
    time = min(time_real,time_energy)

    for t in range(1,time+1): #the real files start 1
        real_name = real_path+'real%04d' % t
        #test#print(real_name)
        real_t = np.loadtxt(real_name)/(-1*au)
        real_im_t = real_t[:,1::2]
        real_re_t = real_t[:,::2]
        energy_t = energy[t-1,:]
        #get Ham_im_
        Ham_im_ = real_re_t
        for i in range(num_orbitals):
            Ham_im_[i,i] = 0
        #get Ham_re_
        Ham_re_ = real_im_t
        for i in range(num_orbitals):
            Ham_re_[i,i] = energy_t[i]
        #test#print(Ham_im_);print(Ham_re_)
        #save Ham_re_ and Ham_im_
        Ham_re_name = real_path+'0_Ham_%s_re' % t
        Ham_im_name = real_path+'0_Ham_%s_im' % t
        np.savetxt(Ham_im_name,Ham_im_,fmt='%.10f')
        np.savetxt(Ham_re_name,Ham_re_,fmt='%.10f')
###################################

###################################
#-task 8
###################################

#####################################################################################################
opts,argvs = getopt.getopt(sys.argv[1:],'j',['task=','nkpt=','band_ini=','band_fin=','time_ini=','time_fin=','a1=','a2=','car=','ispin='])
if opts == []:
    print('*********************************************************************')
    print('--1--Total Energy and Temperature (MD or heating: OSZICAR)         **')
    print('*********************************************************************')
    print('--21-Band Energy (MD: OUTCAR)                                      **')
    print('--22-Average Band Energy (MD: band_energy.dat)                     **')
    print('*********************************************************************')
    print('--3--Standard Deviations of MD trajectory(MD: XDATCAR)             **')
    print('*********************************************************************')
    print('--4--Split XDATCAR into p000X (MD: XDATCAR)                        **')
    print('*********************************************************************')
    print('--51-Distance between A1 and A2 atoms (MD or heating: p000X)       **')
    print('--52-Distance between A1 and A2 atoms (MD or heating: XDATCAR)     **')
    print('*********************************************************************')
    print('--61-Reduce the pristine real files (NAC: real directory:61)       **')
    print('--62-Rename the smaller real files (NAC: newmatrix directory:62)   **')
    print('--63-Reduce the band energy (NAC: band_energy.dat:63)              **')
    print('--64-Average NAC (NAC: smaller newmatrix/real:64)                  **')
    print('--65-Average NAC (NAC: pristine real:65)                           **')
    print('--66-NAC-Time (NAC: pristine real:66)                              **')
    print('*********************************************************************')
    print('--7--Combine Function (NAC: ./res/real000x and ./energy)           **')
    print('*********************************************************************')
    choice = int(input('Which function do you need?----'))
    print('======',choice,'======',choice,'======',choice,'======')
    if choice not in [1,21,22,3,4,5,61,62,63,64,7]:
        print('CHECK INPUT!!!')
    #####################################################################################################

    #####################################################################################################
    ##Total Energy don't need any inputs
    if choice == 1:
        print(' No Input Para!','\n','Need OSZICAR file!','\n','Output TimeTemEn.dat!')
        TotalEne()
    #####################################################################################################
    ##Band Energy need the initial band index and finial band index
    if choice == 21:
        print(' Input nkpt band_ini and band_fin!','\n','Need OUTCAR file!','\n','Output band_energy.dat!')
        ispin = int(input('ispin--'))
        nkpt = int(input('nkpt--'))
        band_ini = int(input('band_ini--'))
        band_fin = int(input('band_fin--'))
        BandEne(task=21,nkpt=nkpt,ispin=ispin,band_ini=band_ini,band_fin=band_fin)
    ##average band energy, print the value from the upper one to lower one
    if choice == 22:
        print(' No Input Para!','\n','Need band_energy.dat!')
        n = 0;average = BandEne(task=22)
        while n < len(average):
            print(round(average[len(average)-n-1],4))
            n+=1
    #####################################################################################################
    ##testing
    ##Standard Deviations of MD trajectory, return a list including standard deviations of each elements.
    if choice == 3:
        print(' Input time_ini and time_fin!','\n','Need XDATCAR file!','\n')
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        save_ =int(input('Save the STD matrix? (Yes:1, No:0)--'))
        StandardDeviations(filename='XDATCAR',save_=save_,time_ini=time_ini,time_fin=time_fin)
    #####################################################################################################
    ##split XDATCAR into p000X
    if choice == 4:
        print(' Input MD_Times!','\n','Need XDATCAR file!','\n','Output pxxxx files!')
        split_xda()
    #####################################################################################################
    ## The distance between A1 and A2 atoms
    if choice == 51:
        print(' Read position from pxxxx files')
        print(' Inputs A1 A2 and pxxxx files num.!','\n','Need pxxxx files!','\n','Output distance.dat!')
        A1 = int(input('The First Atom is:--'))
        A2 = int(input('The Second Atom is:--'))
        time_ini = int(input('time_ini:--'))
        time_fin = int(input('time_fin:--'))
        DistanceTwoAtoms(task=51,A1=A1,A2=A2,Car=1,time_ini=time_ini,time_fin=time_fin)
    if choice == 52:
        print(' Read position from XDATCAR file')
        print(' Inputs A1 A2 and time_ini and time_fin!','\n','Need XDATCAR file!','\n','Output distance.dat!')
        A1 = int(input('The First Atom is:--'))
        A2 = int(input('The Second Atom is:--'))
        time_ini = int(input('time_ini:--'))
        time_fin = int(input('time_fin:--'))
        DistanceTwoAtoms(task=52,A1=A1,A2=A2,Car=0,time_ini=time_ini,time_fin=time_fin)  
    #####################################################################################################
    ##Generally, the NAC matrixs include imaginary part and real part.
    ##They are packed tightly together. Check the shape of real_p.
    ##NOTE: band_ini = VBM-band_min(from NAC calculate)
    ##columns_r=2*(1+band_fin-band_ini)
    if choice == 61:
        print(' Reduce the NAC matrix!')
        print(' Inputs time_fin and time_ini and band_ini and band_fin!','\n','Need realxxxx files!','\n','Output ./newmatrix/realxxxx!')
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        band_ini = int(input('Initial band index is (minus min band in NAC calculate, such as 352-300=52)--'))
        band_fin = int(input('Finial band index is (minus min band in NAC calculate, such as 353-300=53)--'))
        reduce_real(task_=61,band_ini=band_ini,band_fin=band_fin,time_fin=time_fin,time_ini=time_ini)
    if choice == 62:
        print(' Select and Rename the smaller NAC matrix in ./newmatrix/!')
        print(' Inputs time_ini and time_fin!','\n','Need ./newmatrix/realxxxx files!','\n','Output ./500_1500/res/realxxxx!')    
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        reduce_real(task_=62,time_ini=time_ini,time_fin=time_fin)
    if choice == 63:
        print(' Reduce the band_energy.dat!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin(check run_NAC.pbs)!','\n','Need band_energy.dat!','\n','Output reduce_energy!')
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        reduce_real(task_=63,time_ini=time_ini,time_fin=time_fin)
    if choice == 64:
        print(' Average NAC for ./newmatrix/realxxxx!')
        print(' Inputs time_ini and time_fin!','\n','Need ./newmatrix/realxxxx!','\n','Output average_NAC.dat!')
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        reduce_real(task_=64,time_ini=time_ini,time_fin=time_fin)
    if choice == 65:
        print(' Average NAC for pristine realxxxx!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin(check run_NAC.pbs)!','\n','Need pristine realxxxx!','\n','Output average_NAC.dat!')   
        band_ini = int(input('Initial band index is (minus min band in NAC calculate, such as 352-300=52)--'))
        band_fin = int(input('Finial band index is (minus min band in NAC calculate, such as 353-300=53)--'))
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        reduce_real(task_=65,time_ini=time_ini,time_fin=time_fin,band_ini=band_ini,band_fin=band_fin)
    if choice == 66:
        print(' NAC-Time','\n',' Note: run in pristine real files!!!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin!','\n','Need ./newmatrix/realxxxx!','\n','Output NAC_Time.dat!')
        band_ini = int(input('Initial band index is (minus min band in NAC calculate, such as 352-300=52)--'))
        band_fin = int(input('Finial band index is (minus min band in NAC calculate, such as 352-300=52)--'))
        time_ini = int(input('time_ini--'))
        time_fin = int(input('time_fin--'))
        reduce_real(task_=66,time_ini=time_ini,time_fin=time_fin,band_ini=band_ini,band_fin=band_fin)
    #####################################################################################################
    ##combine real and energt to 0_Ham_
    if choice == 7:
        print(' No Inputs!','\n','Need ./energy and ./res/realxxx!','Output 0_Ham_ and energyx!')
        VaspComPyxaidEneNac(real_path='./res/', energy_file='./energy', name_num=4)
    #####################################################################################################
else:
    print(opts)
    for opt,value in opts:
        if opt == '--task':
            task = int(value)
        elif opt == '--nkpt':
            nkpt = int(value)
        elif opt == '--ispin':
            ispin = int(value)
        elif opt == '--band_ini':
            band_ini = int(value)
        elif opt == '--band_fin':
            band_fin = int(value)
        elif opt == '--time_ini':
            time_ini = int(value)
        elif opt == '--time_fin':
            time_fin = int(value)
        elif opt == '--a1':
            a1 = int(value)
        elif opt == '--a2':
            a2 = int(value)
        elif opt == '--car':
            car = int(value)
            
    ##Total Energy don't need any inputs
    if task == 1:
        print(' No Input Para!','\n','Need OSZICAR file!','\n','Output TimeTemEn.dat!')
        TotalEne()
    ##Band Energy need the initial band index and finial band index
    if task == 21:
        print(' nkpt band_ini and band_fin','\n','Need OUTCAR file!','\n','Output band_energy.dat!')
        BandEne(task=21,nkpt=nkpt,band_ini=band_ini,band_fin=band_fin,ispin=ispin)
    if task == 22:
        print(' No Input Para!','\n','Need band_energy.dat!')
        n = 0;average = BandEne(task=22)
        while n < len(average):
            print(round(average[len(average)-n-1],4))
            n+=1
    ##testing
    ##Standard Deviations of MD trajectory, return a list including standard deviations of each elements.
    if task == 3:
        print(' time_ini and time_fin!','\n','Need XDATCAR file!','\n')
        StandardDeviations(filename='XDATCAR',save_=0,time_ini=time_ini,time_fin=time_fin)
    ##split XDATCAR into p000X
    if task == 4:
        print(' Input MD_Times!','\n','Need XDATCAR file!','\n','Output pxxxx files!')
        split_xda(filename='XDATCAR')
    ## The distance between A1 and A2 atoms
    if task == 51:
        print(' Read position from pxxxx files')
        print(' a1 a2 and pxxxx files num.!','\n','Need pxxxx files!','\n','Output distance.dat!')
        DistanceTwoAtoms(task=task,A1=a1,A2=a2,Car=1,time_ini=time_ini,time_fin=time_fin)
    if task == 52:
        print(' Read position from XDATCAR file')
        print(' a1 a2 and time_ini and time_fin!','\n','Need XDATCAR file!','\n','Output distance.dat!')
        DistanceTwoAtoms(task=task,A1=a1,A2=a2,Car=0,time_ini=time_ini,time_fin=time_fin) 
    ##Generally, the NAC matrixs include imaginary part and real part.
    ##They are packed tightly together. Check the shape of real_p.
    ##NOTE: band_ini = VBM-band_min(from NAC calculate)
    ##columns_r=2*(1+band_fin-band_ini)
    if task == 61:
        print(' Reduce the NAC matrix!')
        print(' time_fin and time_ini and band_ini and band_fin!','\n','Need realxxxx files!','\n','Output ./newmatrix/realxxxx!')
        reduce_real(time_fin=time_fin,task_=61,band_ini=band_ini,band_fin=band_fin,time_ini=time_ini)  
    if task == 62:
        print(' Select and Rename the smaller NAC matrix in ./newmatrix/!')
        print(' Inputs time_ini and time_fin!','\n','Need ./newmatrix/realxxxx files!','\n','Output ./500_1500/res/realxxxx!')
        reduce_real(task_=62,time_ini=time_ini,time_fin=time_fin)
    if task == 63:
        print(' Reduce the band_energy.dat!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin(check run_NAC.pbs)!','\n','Need band_energy.dat!','\n','Output reduce_energy!')
        reduce_real(task_=63,time_ini=time_ini,time_fin=time_fin)
    if task == 64:
        print(' Average NAC for ./newmatrix/realxxxx!')
        print(' Inputs time_ini and time_fin!','\n','Need ./newmatrix/realxxxx!','\n','Output average_NAC.dat!')
        reduce_real(task_=64,time_ini=time_ini,time_fin=time_fin)
    if task == 65:
        print(' Average NAC for pristine realxxxx!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin(check run_NAC.pbs)!','\n','Need pristine realxxxx!','\n','Output average_NAC.dat!')   
        reduce_real(task_=65,time_ini=time_ini,time_fin=time_fin,band_ini=band_ini,band_fin=band_fin)
    if task == 66:
        print(' NAC-Time','\n',' Note: run in pristine real files!!!')
        print(' Inputs time_ini and time_fin and band_ini and band_fin!','\n','Need ./newmatrix/realxxxx!','\n','Output NAC_Time.dat!')
        reduce_real(task_=66,time_ini=time_ini,time_fin=time_fin,band_ini=band_ini,band_fin=band_fin)
    ##combine real and energt to 0_Ham_
    if task == 7:
        print(' No Inputs!','\n','Need ./energy and ./res/realxxx!','Output 0_Ham_ and energyx!')
        VaspComPyxaidEneNac(real_path='./res/', energy_file='./energy', name_num=4)
