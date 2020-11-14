# coding: utf-8
#! /usr/bin/python3

'''
    Plot NAC Map
    Input File: NAC matrix and energy file
    Author: yonghao_zhu@163.com
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import seaborn as sns
####################################################################
'''
    Note:
    1. NAC(meV), energy(eV)
    2. Diagonal elements of NAC matrix = 0 meV
    3. Colormap(cmap_, 164)-->https://blog.csdn.net/sinat_32570141/article/details/105356507
    4. deltaE is necessary for BandEne mode, set 0 to use default value.
    5. The band index of NAC file is the same with energy file for BandEne plot mode! 
'''
Plot_NAC_MapFun = {}
Plot_NAC_MapFun['TorF']     = 'T' #T or F
#file sets
Plot_NAC_MapFun['NAC_file'] = 'real0060_300_400' #average NAC file (whole size) 
Plot_NAC_MapFun['ene_file'] = 'band_energy.dat' #band energy file  (whole size)
#parameters set
Plot_NAC_MapFun['band_ini'] = 36 #Include this point (start 0)
Plot_NAC_MapFun['band_fin'] = 42 #Include this point (start 0)
Plot_NAC_MapFun['time_ini'] = 1000 #Include this point (start 0)
Plot_NAC_MapFun['time_fin'] = 1200 #Include this point (start 0)
#show and save set
Plot_NAC_MapFun['show_']    = 'T'
Plot_NAC_MapFun['save_']    = 'F' 
#PlotMode set
Plot_NAC_MapFun['PlotMode'] = 'BandInd' #BandInd or BandEne
Plot_NAC_MapFun['deltaE']   = 0.002 #necessary for BandEne mode 
#color bar sets
Plot_NAC_MapFun['cmap_']    = 'Reds' #164 options
Plot_NAC_MapFun['CBar_max'] = 2 #set the value 
Plot_NAC_MapFun['CBar_min'] = 0 #0 is OK for NAC!
Plot_NAC_MapFun['CBar_lab'] = 5 #the number of lables of color bar
####################################################################
def Plot_NAC_Map(NAC_file='real0060_300_400', ene_file='band_energy.dat', band_ini=1,band_fin=2, time_ini=1,time_fin=2,
    show_='T',save_='T',plotmode='BandInd',labelsize=15,cmap_='Reds',
    deltaE=0,CBar_max=1,CBar_min=0,CBar_lab=5):
    if os.path.isfile(NAC_file):
        nac_pri = np.loadtxt(NAC_file)
        print('--------------Matrix Shape--------------')
        print('NAC matrix shape:',nac_pri.shape)

        ##Get Data
        #reduce nac matrix, remove imaginary part
        nac_matrix = 1000*abs(nac_pri[band_ini:band_fin+1,2*band_ini:2*band_fin+2][:,::2])
        #not reduce
        #nac_matrix = 1000*nac_pri
        #set diagonal elements of NAC matrix to 0 meV

        row,col = np.diag_indices_from(nac_matrix)
        nac_matrix[row,col] = 0
        #test#print(nac_matrix.shape)

        #reduce energy matrix
        if os.path.isfile(ene_file):
            ene_pri = np.loadtxt(ene_file)
            print('Ene matrix shape:',ene_pri.shape)
            print('----------------------------------------')
            ene = np.mean(ene_pri[time_ini:time_fin+1,band_ini:band_fin+1],axis=0)
        #test#print(ene)
        else:
            print('No Energy File!')
        ##plot NAC map
        #font-->Times New Roman
        plt.rc('font',family='Times New Roman')
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        #x,y-->band index
        if plotmode == 'BandInd':
            fig,ax = plt.subplots()
            htmap=sns.heatmap(nac_matrix,square=True,cmap=cmap_,vmax=CBar_max,vmin=CBar_min,
                linewidths=0.2,cbar=False) #linecolor='reds'
            #color bar set
            cbar=htmap.figure.colorbar(htmap.collections[0])

            cbar_lable=['%.2f' % (CBar_min+i*(CBar_max-CBar_min)/(CBar_lab-1)) for i in range(CBar_lab)]
            cbar.set_ticks([i*(CBar_max-CBar_min)/(CBar_lab-1) for i in range(CBar_lab)])

            cbar.set_ticklabels(cbar_lable)
            cbar.ax.tick_params(labelsize=15)
            cbar.set_label('meV',fontsize=20)
            #plot set
            ax.invert_yaxis()
            plt.xlabel('Band Index',size=20); plt.ylabel('Band Index',size=20) 
            plt.xticks([i+0.5 for i in range(nac_matrix.shape[0])],[i+1 for i in range(nac_matrix.shape[0])],size=labelsize)
            plt.yticks([i+0.5 for i in range(nac_matrix.shape[0])],[i+1 for i in range(nac_matrix.shape[0])],size=labelsize)
            fig.subplots_adjust(top=0.9, bottom=0.2, left=0.1, right=0.9)
            plt.tick_params(bottom=False,left=False)
            if save_ == 'T':
                #plt.savefig('NAC-Map-BI.svg',format='svg')
                plt.savefig("NAC-Map-BI.png",dpi=500,bbox_inches = 'tight')
        #x,y-->band energy
        if plotmode == 'BandEne':
            #insert points in NAC matrix
            if deltaE != 0:
                insert_num = [int((ene[i+1]-ene[i])/deltaE)+1 for i in range(len(ene)-1)]
            if deltaE == 0:
                deltaE = min([ene[i+1]-ene[i] for i in range(len(ene)-1)])
                print('Default deltaE=',deltaE)
                insert_num = [int((ene[i+1]-ene[i])/deltaE)+1 for i in range(len(ene)-1)]
            print('insert_num-->',insert_num)
            print('old_nac_matrix.shape-->',nac_matrix.shape)
            size_new_nac = sum(insert_num)+nac_matrix.shape[0]
            print('new_nac_matrix.shape--> (',size_new_nac,',',size_new_nac,')')
            print('----------------------------------------')

            #insert row
            '''
            0 2-->0 *0.5* *1.0* *1.5* 2
            2 0-->2 *1.5* *1.0* *0.5* 0
            '''
            new_nac_matrix = []
            for row in range(nac_matrix.shape[0]): 
                x0 = []
                for col in range(nac_matrix.shape[1]-1):
                    x1 = [nac_matrix[row,col]
                        +i*(nac_matrix[row,col+1]-nac_matrix[row,col])/(insert_num[col]+1)
                        for i in range(insert_num[col]+1)]
                    for row_new in x1:
                        x0.append(row_new)
                x0.append(nac_matrix[row,-1])
                new_nac_matrix.append(x0)
            new_nac_matrix = np.array(new_nac_matrix)

            #build new_nac.shape
            new_nac = np.zeros([size_new_nac,size_new_nac])
            row_0 = list(new_nac_matrix[0])
            new_nac[0,:] = row_0; new_nac[:,0] = row_0
            row_1 = list(new_nac_matrix[-1])
            new_nac[-1,:] = row_1; new_nac[:,-1] = row_1
            insert_num_sum = [sum(insert_num[:i+1]) for i in range(len(insert_num))]
            #test#print(insert_num_sum)#[2, 7, 9, 15, 21, 29]
            #test#print(new_nac[-80:-1,-1])
            for i in range(1,len(insert_num_sum)):
                new_nac[insert_num_sum[i-1]+i,:] = new_nac_matrix[i]
                new_nac[:,insert_num_sum[i-1]+i] = new_nac_matrix[i]
            m = new_nac
            #test#print(new_nac[-80:-1,-1])
            #write new_nac-->step1
            insert_num_sum = [0]+insert_num_sum
            temp = []
            for col in range(len(insert_num_sum)-1):
                temp_1 = insert_num_sum[col]+col
                temp_2 = insert_num_sum[col+1]+col+1
                temp.append([temp_1,temp_2])
            #test#print(temp)#[[0, 3], [3, 9], [9, 12], [12, 19], [19, 26], [26, 35]]
            for row in range(new_nac.shape[0]):
                if list(new_nac[row,:]).count(0) > 1: #the number of 0 element
                    for ii in temp:
                        if row < ii[1]:
                            mat0 = list(new_nac[row,ii[0]:ii[1]+1])
                            if row > ii[0]:
                                mat1 = [0+j*(mat0[ii[1]-ii[0]]-0)/(ii[1]-row) for j in range(ii[1]-row+1)]
                                new_nac[row,row:ii[1]+1] = mat1
                            if row <= ii[0]:
                                mat1 = [mat0[0]+j*(mat0[-1]-mat0[0])/(len(mat0)-1) for j in range(len(mat0))]
                                new_nac[row,ii[0]:ii[1]+1] = mat1
            #write new_nac-->step2
            for i in range(1,new_nac.shape[0]):
                for j in range(1,new_nac.shape[1]):
                    if j > i:
                        new_nac[j,i] = new_nac[i,j] 
            #test#print(new_nac.shape)
            fig,ax = plt.subplots()
            htmap=sns.heatmap(new_nac,square=True,cmap=cmap_,vmax=CBar_max,vmin=CBar_min,
                cbar=False)
            #color bar set
            cbar=htmap.figure.colorbar(htmap.collections[0])
            cbar_lable=['%.2f' % (CBar_min+i*(CBar_max-CBar_min)/(CBar_lab-1)) for i in range(CBar_lab)]
            cbar.set_ticks([CBar_min+i*(CBar_max-CBar_min)/(CBar_lab-1) for i in range(CBar_lab)])
            cbar.set_ticklabels(cbar_lable)
            cbar.ax.tick_params(labelsize=15)
            cbar.set_label('meV',fontsize=20)
            ax.invert_yaxis()
            #plot set
            plt.xlabel('Band Energy',size=20); plt.ylabel('Band Energy',size=20) 
            plt.xticks([i*int(new_nac.shape[0]/4) for i in range(5)],['%.2f' % (min(ene)+i*((max(ene)-min(ene))/4)) for i in range(5)],size=labelsize)
            plt.yticks([i*int(new_nac.shape[0]/4) for i in range(5)],['%.2f' % (min(ene)+i*((max(ene)-min(ene))/4)) for i in range(5)],size=labelsize)
            plt.xticks(rotation=45); plt.yticks(rotation=45)
            fig.subplots_adjust(top=0.9, bottom=0.22, left=0.1, right=0.9)
            plt.tick_params(bottom=False,left=False)

            if save_ == 'T':
                #plt.savefig('NAC-Map-BE.svg',format='svg') 
                plt.savefig("NAC-Map-BE.png",dpi=500,bbox_inches = 'tight')

        if show_ == 'T':
            plt.show()
    else:
        print('No Files!!!')

if Plot_NAC_MapFun['TorF'] == 'T':
    Plot_NAC_Map(NAC_file=Plot_NAC_MapFun['NAC_file'],
        ene_file=Plot_NAC_MapFun['ene_file'],
        band_ini=Plot_NAC_MapFun['band_ini'],
        band_fin=Plot_NAC_MapFun['band_fin'],
        time_ini=Plot_NAC_MapFun['time_ini'],
        time_fin=Plot_NAC_MapFun['time_fin'],
        show_=Plot_NAC_MapFun['show_'],save_=Plot_NAC_MapFun['save_'],cmap_=Plot_NAC_MapFun['cmap_'],
        plotmode=Plot_NAC_MapFun['PlotMode'],
        deltaE=Plot_NAC_MapFun['deltaE'],
        CBar_max=Plot_NAC_MapFun['CBar_max'],
        CBar_min=Plot_NAC_MapFun['CBar_min'],
        CBar_lab=Plot_NAC_MapFun['CBar_lab'])