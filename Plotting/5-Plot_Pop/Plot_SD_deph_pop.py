# coding: utf-8
#！/usr/bin/python3

'''
    Plot Spectral-Density and fit Dephsing-Time with pyxaid format for test tasks
    author: Yonghao Zhu (yonghao_zhu@163.com)
'''

################################################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
################################################

##################FitDephasing##################

FitDephasing_fun                   = {}
#FitGauss->y=a*np.exp(-0.5*(x/b)**2)
FitDephasing_fun['TorF']           = True # True or False
FitDephasing_fun['mode']           = 'pyxaid' # pyxaid or others
# pyxiad
FitDephasing_fun['icond']          = '5'
FitDephasing_fun['pair']           = '0_1'
FitDephasing_fun['Dt']             = False # True or False
# others
FitDephasing_fun['file']           = 'auto_dephasing.dat' 

FitDephasing_fun['Time_cutoff']    = 200
FitDephasing_fun['show_']          = True # True or False
#####################PlotSD######################
PlotSD_fun                 = {}
PlotSD_fun['TorF']         = True # True or False
PlotSD_fun['mode']         = 'pyxaid' # pyxaid or others
# pyxaid
PlotSD_fun['icond']        = '5'
PlotSD_fun['pair']         = '0_1'
#others
PlotSD_fun['file_']        = 'sd.dat'
PlotSD_fun['Freq_cutoff']  = 1000
PlotSD_fun['show_']        = True # True or False
#####################PlotPOP#####################
PlotPop_fun                 = {}
'''
    only two states
'''
PlotPop_fun['TorF']         = True # True or False
PlotPop_fun['PlotType']     = 'libra' # v8 (v8+hefei_namd), libra, pyxaid,
# v8        : out.avg
# libra     : _out.txt
# hefei_namd: Ave_SHPROP.dat
# pyxaid    : jj-pop.dat
PlotPop_fun['file_name']    = '0_out.txt'
#-------------------------------------------------
# libra plot
PlotPop_fun['dt_libra']     = 41 # in a.u. of time
#-------------------------------------------------
# pyxaid, v8 and hefei: -1 is all, s0->0, s1->1
# libra: 
    # nstates = 2 => n_cols = 2*3 + 5 = 10
    #        0    1      2     3    4      5      6     7     8    9   10
    # res =  t  E0, P_SE0, P_SH0, E1, P_SE1, P_SH1, E_SE, E_SH, P_SE, P_SH
PlotPop_fun['Sele_states']  = 3 

PlotPop_fun['Time_cutoff']  = 5000 #fs
PlotPop_fun['show_']        = True # True or False
#FitLine->y=a*x+b
#FitExponential->y=y0+a1*exp(-1*x/t1)
PlotPop_fun['Fit_']         = 'Line' #Line or Exponential function
PlotPop_fun['Stage']        = 1 #1 or 2, for Exponential Fitting
#################################################

def FitGauss(x,a,b):
    x = np.array(x)
    return a*np.exp(-0.5*(x/b)**2)

def FitLine(x,a,b):
    x = np.array(x)
    return a*x+b

def FitExponential_1(x,y0,a1,t1):
    x = np.array(x)
    return y0+a1*np.exp(-1*x/t1)

def FitExponential_0(x,t):
    x = np.array(x)
    return 1-np.exp(-1*x/t)

#################################################
def FitDephasing(file_='auto_dephasing.dat', mode='pyxaid', Dt=False, Time_cutoff=200, show_=False):

    if mode == 'pyxaid':

        dephasing_data = np.loadtxt(file_,comments='Time')
        x_time = dephasing_data[:,0][:Time_cutoff]

        if Dt:
            y_dephasing = dephasing_data[:,1][:Time_cutoff]
        else:
            y_dephasing = dephasing_data[:,2][:Time_cutoff]

    if mode == 'others':
        dephasing_data = np.loadtxt(file_)
        x_time = dephasing_data[:,0][:Time_cutoff]

        y_dephasing = dephasing_data[:,2][:Time_cutoff]

    para,pcov = curve_fit(FitGauss,x_time,y_dephasing)

    print('Dephasing Time = ', round(para[1],4),'fs')

    if show_:
        plt.rc('font',family='Times New Roman')
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        a = para[0];b = para[1] 
        y_fit = FitGauss(x_time,a,b)
        plt.scatter(x_time,y_dephasing)
        plt.plot(x_time,y_fit,linewidth=1)

    if show_:
        plt.xlim(0,Time_cutoff);plt.ylim([0,1])
        plt.ylabel("Dephasing",size=20);plt.xlabel("Time (fs)",size=20)
        plt.show()

if FitDephasing_fun['TorF']:

    if FitDephasing_fun['mode'] == 'pyxaid':

        file_name = 'icond' + FitDephasing_fun['icond'] + 'pair' + FitDephasing_fun['pair'] + 'Dephasing_function.txt'

        if os.path.isfile(file_name):
            FitDephasing(file_ = file_name, mode='pyxaid', Dt=FitDephasing_fun['Dt'], 
                         Time_cutoff=FitDephasing_fun['Time_cutoff'], show_=FitDephasing_fun['show_'])
        else:
            print('No %s!' % file_name)
            exit()

    if FitDephasing_fun['mode'] == 'others':
        if os.path.isfile(FitDephasing_fun['file']):
            FitDephasing(file_ = FitDephasing_fun['file'], mode='others', 
                         Time_cutoff=FitDephasing_fun['Time_cutoff'], show_=FitDephasing_fun['show_'])
        else:
            print('No %s!' % FitDephasing_fun['file'])
            exit()

#################################################
def PlotSD(mode='pyxaid', file_='sd.dat',
           Freq_cutoff=1000, show_=False):
    if mode == 'pyxaid':
        SD_data = np.loadtxt(file_,dtype=str)[:,1::2][:,1:3].astype(np.float)
        x_freq = SD_data[:,0]
        y_inten = SD_data[:,1]
        
    if mode == 'others':
        SD_data = np.loadtxt(file_)
        x_freq = SD_data[:,0]
        y_inten = SD_data[:,1]

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in' 
    plt.plot(x_freq,y_inten,linewidth=2)
    plt.xlim(0,Freq_cutoff)   

    plt.ylabel("Spectral Density",size=20);plt.xlabel("Frequency ($cm^{-1}$)",size=20)
    plt.ylim([0,max(y_inten)])
    
    if show_:
        plt.show()

if PlotSD_fun['TorF']:

    if PlotSD_fun['mode'] == 'pyxaid':
        file_name = 'icond' + PlotSD_fun['icond'] + 'pair' + PlotSD_fun['pair'] + 'Spectral_density.txt'

        if os.path.isfile(file_name):
            PlotSD(mode='pyxaid', file_=file_name,
                   Freq_cutoff=PlotSD_fun['Freq_cutoff'], show_=PlotSD_fun['show_'])
        else:
            print('No %s!' % file_name)
            exit()

    if PlotSD_fun['mode'] == 'others':

        if os.path.isfile(PlotSD_fun['file_']):
            PlotSD(mode='others', file_=PlotSD_fun['file_'],
                   Freq_cutoff=PlotSD_fun['Freq_cutoff'], show_=PlotSD_fun['show_'])            
        else:
            print('No %s!' % PlotSD_fun['file_'])
            exit()

#################################################
def PlotPop(Time_cutoff=1000, PlotType='pyaxid', show_=False, fit_='Line', Sele_states=-1, Stage=1, dt_libra=41):

    if PlotType == 'pyxaid': #rates and population fitting 
        pop_data = np.loadtxt(PlotPop_fun['file_name'])[0:Time_cutoff,:]
        num_states = pop_data.shape[1]-1
        y_pop = pop_data[:,1:]

    if  PlotType == 'v8': #rates fitting
        pop_data = np.loadtxt(PlotPop_fun['file_name'],comments=['t'])[0:Time_cutoff,:]
        num_states = pop_data.shape[1]-4
        y_pop = pop_data

    if PlotType == 'libra': #rates and population fitting 
        pop_data = np.loadtxt(PlotPop_fun['file_name'])[0:Time_cutoff,:]
        num_states = int((pop_data.shape[1]-5)/3)
        y_pop = pop_data

    if PlotType in ['v8', 'pyxaid']:
        x_time = pop_data[:,0]
    if PlotType in ['libra']:
        x_time = pop_data[:,0]/dt_libra
    #test#print(x_time.shape); print(y_pop.shape)

    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    if Sele_states == -1:#All states->plot
        y_max,y_min = [],[]
        for i in range(num_states):
            y_max.append(max(y_pop[:,i])); y_min.append(min(y_pop[:,i]))
            if PlotType == 'pyxaid':#pyxaid
                plt.scatter(x_time/1000,y_pop[:,i],s=10)
            if PlotType == 'v8':#v8
                plt.scatter(x_time,y_pop[:,i],s=10)
        if show_:
            if PlotType == 'pyxaid':#pyxaid
                plt.xlim(0,Time_cutoff/1000); plt.xticks(size=15)
                plt.ylabel("Population",size=20);plt.xlabel("Time (ps)",size=20)
            if PlotType == 'v8':#v8
                plt.xlim(0,Time_cutoff); plt.xticks(size=15)
                plt.ylabel("Population",size=20);plt.xlabel("Time (fs)",size=20)
            plt.ylim(min(y_min),max(y_max)); plt.yticks(size=15)
            plt.show()

    if Sele_states != -1:#Plot and Fits
        
        if PlotType == 'pyxaid':#pyxaid
            y_pop_ = y_pop[:,Sele_states]
            plt.scatter(x_time/1000,y_pop_,s=10)
            plt.xlim(0,Time_cutoff/1000); plt.xticks(size=15)
            plt.ylabel("Population",size=20);plt.xlabel("Time (ps)",size=20)

        if PlotType == 'v8':#v8
            if Sele_states == 0:
                y_pop_ = y_pop[:,-2]
            if Sele_states == 1:
                y_pop_ = y_pop[:,-1]
            plt.scatter(x_time,y_pop_,s=10)
            plt.xlim(0,Time_cutoff); plt.xticks(size=15)
            plt.ylabel("Population",size=20);plt.xlabel("Time (fs)",size=20)

        if PlotType == 'libra':#libra
            y_pop_ = y_pop[:,Sele_states]
            plt.scatter(x_time,y_pop_,s=10)
            plt.xlim(0,Time_cutoff); plt.xticks(size=15)
            plt.ylabel("Population",size=20);plt.xlabel("Time (fs)",size=20)

        plt.ylim(min(y_pop_),max(y_pop_)); plt.yticks(size=15)

        #fitting
        if fit_ == 'Line':
            para,pcov = curve_fit(FitLine,x_time,y_pop_)
            print('Rates (s%s, fs^-1) = ' % Sele_states,para[0],'fs^-1')
            print('Population Time (s%s, ps) = ' % Sele_states,round(0.001/para[0],8),'ps')
            y_fit = FitLine(x_time,para[0],para[1])
            if PlotType == 'pyxaid':#pyxaid
                plt.plot(x_time/1000,y_fit,linewidth=2,linestyle='--',color='red')
            if PlotType in ['v8', 'libra']:#v8 and libra
                plt.plot(x_time,y_fit,linewidth=2,linestyle='--',color='red')

        if fit_ == 'Exponential':
            if Stage == 1:
                if Sele_states != 0:
                    para,pcov = curve_fit(FitExponential_1,x_time,y_pop_)
                    print('Population Time (s%s, ps) = ' % Sele_states,round(para[2]/1000,8),'ps')
                    y_fit = FitExponential_1(x_time,para[0],para[1],para[2])
                if Sele_states == 0:
                    para,pcov = curve_fit(FitExponential_0,x_time,y_pop_)
                    print('Population Time (s%s, ps) = ' % Sele_states,round(para[0]/1000,8),'ps')
                    y_fit = FitExponential_0(x_time,para[0])
                #test#print(para[0]); print(para[1]); print(para[2])
                if PlotType == 'pyxaid':#pyxaid
                    plt.plot(x_time/1000,y_fit,linewidth=2,linestyle='--',color='red')
                if PlotType in ['v8', 'libra']:#v8 and libra
                    plt.plot(x_time,y_fit,linewidth=2,linestyle='--',color='red')

            if Stage == 2:
                y_pop_max = np.max(y_pop_)
                max_index = np.where(y_pop_==y_pop_max)[0][-1]
                #first stage
                y_pop_1 = y_pop_[:max_index]
                x_time_1 = x_time[:max_index]
                para_1,pcov_1 = curve_fit(FitExponential_0,x_time_1,y_pop_1)
                print('Population Time (s%s, first stage, ps) = ' % Sele_states,round(para_1[2]/1000,8),'ps')
                y_fit_1 = FitExponential_0(x_time_1,para_1[0],para_1[1],para_1[2])
                if PlotType == 'pyxaid':#pyxaid
                    plt.plot(x_time_1/1000,y_fit_1,linewidth=2,linestyle='--',color='red')
                if PlotType == 'v8':#v8
                    plt.plot(x_time_1,y_fit_1,linewidth=2,linestyle='--',color='red')

                #sencond stage
                y_pop_2 = y_pop_[max_index:]
                x_time_2 = x_time[max_index:]-max_index
                para_2,pcov_2 = curve_fit(FitExponential_1,x_time_2,y_pop_2)
                print('Population Time (s%s, second stage, ps) = ' % Sele_states,round(para_2[2]/1000,8),'ps')
                y_fit_2 = FitExponential_1(x_time_2,para_2[0],para_2[1],para_2[2])
                if PlotType == 'pyxaid':#pyxaid
                    plt.plot(x_time_2/1000+max_index/1000,y_fit_2,linewidth=2,linestyle='--',color='black')
                if PlotType == 'v8':#v8
                    plt.plot(x_time_2+max_index,y_fit_2,linewidth=2,linestyle='--',color='black')

        if show_:
            plt.show() 


if PlotPop_fun['TorF']:
    if os.path.isfile(PlotPop_fun['file_name']):
        PlotPop(Time_cutoff=PlotPop_fun['Time_cutoff'], PlotType=PlotPop_fun['PlotType'],
            show_=PlotPop_fun['show_'], fit_=PlotPop_fun['Fit_'],
            Sele_states=PlotPop_fun['Sele_states'],
            Stage=PlotPop_fun['Stage'], dt_libra=PlotPop_fun['dt_libra'])
    else:
        print('No %s!' % PlotPop_fun['file_name'])

#################################################