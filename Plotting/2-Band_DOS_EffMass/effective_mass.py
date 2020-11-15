#! /usr/bin/python3
'''
    This script can read band.dat for effective mass.
    The file, band.dat obtained from vaspkit is necessary.
    NOTE:
    (1) the kpt is different between vaspkit and p4vasp format, vaspkit = p4*2pi
    (2) 1eV = 0.036749 Hartree;1Bohr = 0.529177 Ang
    efftive mass = 0.5/(c*0.036749/0.529177/0.529177)
    (3) Check the vaspkit version if error.
'''

import os
import numpy as np
from scipy.optimize import leastsq
import linecache
import matplotlib as mpl
import matplotlib.pyplot as plt
############################################
EffMassFun = {}
#input files: band.dat BAND_GAP KLABELS
#The number of VBM/CBM in a high kpoints path is up to 1! 
EffMassFun['TorF']     = 'T'
#0: normal->read bands from band.dat
#1: unormal->read bands from a band level for more than one VBM/VBM point.
EffMassFun['Type']     = 0 #0 or 1
EffMassFun['CorV']     = 'V' #C->Conduction Band V->Valence Band
EffMassFun['iband']    = 10 #start 1
EffMassFun['Fit_nkpt'] = 16
#Plot
EffMassFun['show_']    = 'T'
############################################

def bandKpt(iband=11):
    '''
    file=Band.dat, BAND_GAP, KLABELS
    	#K-Path   Energy-Level
    	# NKPTS & NBANDS:  80 24
    iband = 11 #start 1
    '''
    line_2 = linecache.getline('band.dat',2).split()
    hkpt = [i.split()[1] for i in linecache.getlines('KLABELS') if len(i.split()) == 2]
    nkpts,nbands = int(line_2[-2]),int(line_2[-1])
    band = np.loadtxt('band.dat',comments=['#']).reshape(nbands,nkpts,2)
    #test#print(band.shape)#print(band[0])
    num_vb = int(linecache.getline('BAND_GAP',7).split()[4]); num_cb = nbands-num_vb
    num_hkpt = len(hkpt)-1
    ene = band[iband-1,:,1].reshape(num_hkpt,-1)
    kpt = band[iband-1,:,0].reshape(num_hkpt,-1)
    return ene,kpt,band

class fit:
    def __init__(self):
        self.para = [];self.r2 = [] 
    def func(self,params,x):
        a,b,c = params
        return a*x*x+b*x+c
    def error(self,params,x,y):
        return self.func(params,x)-y
    def SlovePara(self,m):
        X = m[0];Y = m[1];p0 = [10,10,10]
        Para = leastsq(self.error,p0,args=(X,Y))
        self.para = Para[0];f = np.array([self.para[0]*i*i+self.para[1]*i+self.para[2] for i in X])
        self.r2 = 1-(np.sum(np.square(Y-f)))/(np.sum(np.square(Y-np.mean(Y))))
##test fit############################################################################################
# m =  np.array([[0.000,0.029,0.057,0.086,0.114,0.143],[0.1599,0.1799,0.2362,0.3202,0.4220,0.5341]])##
# b = fit()                                                                                         ##
# b.SlovePara(m)                                                                                    ##
# print(b.r2)                                                                                       ##
######################################################################################################

def calculate_fit(ene,kpt,CorV='C',Fit_nkpt=6):
    '''
    CorV=C or V

    '''
    if CorV == 'C':##CB
        min_path = np.argwhere(ene == np.min(ene))[:,0]
        
    if CorV == 'V':##VB
        min_path = np.argwhere(ene == np.max(ene))[:,0]

    #test#print(min_path)#print(len(ene))

    if len(min_path) >= 2:
        return []
        print('ERROR')          

    if len(min_path) == 1:
        if min_path[0] == len(ene)-1:#right
            X = [kpt[len(ene)-1-i] for i in range(Fit_nkpt)]; Y = [ene[len(ene)-1-i] for i in range(Fit_nkpt)]
            fitting = fit();fitting.SlovePara(np.vstack((X,Y)));a = fitting.para[0];r2 = fitting.r2
            return a,r2,X,Y

        if min_path[0] == 0:#left
            X = [kpt[i] for i in range(Fit_nkpt)]; Y = [ene[i] for i in range(Fit_nkpt)]
            fitting = fit();fitting.SlovePara(np.vstack((X,Y)));a = fitting.para[0];r2 = fitting.r2
            return a,r2,X,Y

        if min_path[0] > 0 and min_path < len(ene)-1:#middle
            half = int(Fit_nkpt/2)
            if half-1 < min_path[0] and min_path[0] < len(ene)-half+1:
                X = [kpt[i+min_path[0]] for i in range(-1*half,half+1)]
                Y = [ene[i+min_path[0]] for i in range(-1*half,half+1)]

            if min_path[0] <= half-1:
                X = [kpt[i+min_path[0]] for i in range(-1*min_path[0],half)]
                Y = [ene[i+min_path[0]] for i in range(-1*min_path[0],half)]

            if min_path[0] >= len(ene)-half+1:
                X = [kpt[i+min_path[0]] for i in range(-1*half,len(ene)-min_path[0])]
                Y = [ene[i+min_path[0]] for i in range(-1*half,len(ene)-min_path[0])]

            fitting = fit();fitting.SlovePara(np.vstack((X,Y)));a = fitting.para[0];r2 = fitting.r2
            return a,r2,X,Y

def EffMass(iband=11,CorV='C',Fit_nkpt=6,Type=0,
    show_='T'):
    '''
    file=Band.dat
    	#K-Path   Energy-Level
    	# NKPTS & NBANDS:  80 24
    CorV=C or V
    '''
    if Type == 0:
        ene,kpt,band = bandKpt(iband=iband)
    else:
        print('Sorry, this funciotn is on the way!')
    #test#print(len(ene),len(kpt))
    print('**********OUTPUT**********')
    plot_data = []
    for i in range(len(ene)):
        a,r2,X,Y = calculate_fit(ene=ene[i],kpt=kpt[i],Fit_nkpt=Fit_nkpt,CorV=CorV) 
        plot_data.append([a,X,Y])
        '''
        1eV = 0.036749 Hartree;1Bohr = 0.529177 Ang
        '''
        if CorV == 'C':
            ele_eff = round(0.5/(a*0.036749/0.5291772108/0.5291772108),4)
            print('CBM:',min(Y),'eV')
            print('the fitting parameter, a:','\n',round(a,4))
            print('the fitting parameter, r2:','\n',round(r2,4))
            print('the effective mass of electron:','\n',ele_eff,'me')
            print('------------------------')
        if CorV == 'V':
            print('VBM:',max(Y),'eV')
            print('the fitting parameter, a:','\n',round(a,4))
            print('the fitting parameter, r2:','\n',round(r2,4))
            hole_eff = round(0.5/(a*0.036749/0.5291772108/0.5291772108),4)
            print('the effective mass of electron:','\n',hole_eff,'me')
            print('------------------------')
    print('**********END**********')
    #plot
    if show_ == 'T':
        plt.rc('font',family='Times New Roman')
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        
        for i in range(len(plot_data)): 
            x = plot_data[i][1]
            y = plot_data[i][2]
            band_line = plt.scatter(x,y,linewidth=5)
        plt.ylabel("Energy (eV)",size=20)
        plt.xlabel("kpoints",size=20)
        plt.show()

if EffMassFun['TorF'] == 'T':
	if os.path.isfile('band.dat') and os.path.isfile('BAND_GAP') and os.path.isfile('KLABELS'):	
		EffMass(Fit_nkpt=EffMassFun['Fit_nkpt'],iband=EffMassFun['iband'],Type=EffMassFun['Type'],CorV=EffMassFun['CorV'],
            show_=EffMassFun['show_'])
	else:
		print('Error files')
	