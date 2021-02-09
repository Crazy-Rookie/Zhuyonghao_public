#! /usr/bin/python3
# -*- conding=UTF-8 -*-

'''
    1. print help
    2. author: yonghao_zhu@163.com
'''

def printhelp(task=-1):
    '''
    task = -1: all 
    task = [1, 2, 3, 4, 5, 6, 7]
    '''

    if task in [1, 2, 3, 4, 5, 6, 7]:
        print('Check mdpy help for task %s.' % task)
    elif task == -1:
        print('Check mdpy help for all tasks!')
    else:
        print('Check input!!!')
        exit()

    print('-------------------------------------------------------------')
    print('Hello, I was wirtten for MD preprocessing and postprocessing!')
    print('Author: yonghao_zhu@163.com')
    print('-------------------------------------------------------------')

    if task == 1 or task == -1:
        print('************************************************************')
        print('task --1:         Time-Total_Energy and Temperature')
        print('Input files:      OSZICAR (MD)')
        print('Input parameters: none')
        print('Output files:     TimeTemEn.dat-->Time Total_Energy Temperature')
        print('Usage:            mdpy --task 1')
        print('************************************************************')

    if task == 2 or task == -1:
        print('************************************************************')
        print('task --2*:          Time-BandEnergy')
        print('  task --21:        save Time-BandEnergy')
        print('  Input files:      OUTCAR (MD)')
        print('  Input parameters: nspin nkpt band_ini band_fin')
        print('  Output files:     band_energy.dat--> band1 band2 band3 ...')
        print('  Usage:            mdpy --task 21 --nspin 1 --nkpt 1 --band_ini 31 --band_fin 32')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --22:        average BandEnergy')
        print('  Input files:      band_energy.dat')
        print('  Input parameters: none')
        print('  Output files:     none')
        print('  Usage:            mdpy --task 22')
        print('************************************************************')

    if task == 3 or task == -1:
        print('************************************************************')
        print('task --3:         Standard Deviations')
        print('Input files:      XDATCAR (MD)')
        print('Input parameters: time_ini time_fin')
        print('Ouput files:      STD_atoms.dat')
        print('Usage:            mdpy --task 3 --time_ini 1 --time_fin 1000')
        print('************************************************************')

    if task == 4 or task == -1:
        print('************************************************************')
        print('task --4:         split XDATCAR')
        print('Input files:      XDATCAR (MD)')
        print('Input parameters: MD_times')
        print('Ouput files:      ./pfiles/pxxxx')
        print('Usage:            mdpy --task 4')
        print('************************************************************')

    if task == 5 or task == -1:
        print('************************************************************')
        print('task --5*:          Distance between Two Atoms')
        print('  task --51:        From pfiles')
        print('  Input files:      ./pflies/pxxxx')
        print('  Input parameters: A1 A2 time_ini time_fin')
        print('  Output files:     none')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --52:        From XDATCAR')
        print('  Input files:      XDATCAR')
        print('  Input parameters: A1 A2 time_ini time_fin')
        print('  Output files:     none')
        print('Usage: mdpy --task 51/52 --a1 10 --a2 11 --time_ini 1 --time_fin 1000')
        print('************************************************************')

    if task == 6 or task == -1:
        print('************************************************************')
        print('task --6*:          NACs')
        print('  task --61:        Reduce the pristine real files')
        print('  Input files:      ./real****')
        print('  Input parameters: nac_scale(1 or 0.33) band_ini band_fin time_ini time_fin')
        print('  Ouput files:      ./newmatrix/real****')
        print('  NOTE:             real_old * nac_scale (1 or 0.33) = real_new')
        print('  Usage:            mdpy --task 61 --nac_scale 1 --band_ini 1 --band_fin 2 --time_ini 1 --time_fin 1000')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --62:        Rename the reduced real files')
        print('  Input files:      ./real****')
        print('  Input parameters: time_ini time_fin')
        print('  Ouput files:      ./time_ini_time_fin/res/real****')
        print('  Usage:            mdpy --task 62 --time_ini 1 --time_fin 1000')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --63:        Reduce the band energy')
        print('  Input files:      band_energy.dat')
        print('  Input parameters: time_ini time_fin')
        print('  Ouput files:      reduce_energy')
        print('  Usage:            mdpy --task 63 --time_ini 1 --time_fin 1000')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --64:        Average NAC from reduced real files')
        print('  Input files:      ./real****')
        print('  Input parameters: time_ini time_fin')
        print('  Ouput files:      average.dat')
        print('  Usage:            mdpy --task 64 --time_ini 1 --time_fin 1000')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --65:        Average NAC from pristine real files')
        print('  Input fiels:     ./real****')
        print('  Input parameters: band_ini band_fin time_ini time_fin')
        print('  Ouput files:      average.dat')
        print('  Usage:            mdpy --task 65 --band_ini 10 --band_fin 11 --time_ini 1 --time_fin 1000')
        print('  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('  task --66:        NAC-time')
        print('  Input files:     ./real****')
        print('  Input parameters: band_ini(start 0) band_fin(start 0) time_ini time_fin')
        print('  Ouput files:      NAC_Time.dat')
        print('  Usage:            mdpy --task 66 --band_ini 10 --band_fin 11 --time_ini 1 --time_fin 1000')
        print('************************************************************')

    if task == 7 or task == -1:
        print('************************************************************')
        print('task --7:         combine the energy and nac files to Ham_*_im and Ham_*_re files')
        print('Input files:      ./res/real**** + ./energy')
        print('Input parameters: none')
        print('Ouput files:      ./res/Ham_*_im/res')
        print('Usage:            mdpy --task 7')
        print('***********************************************************')




def main():
    printhelp(task = -1)

if __name__ == "__main__":
    main()