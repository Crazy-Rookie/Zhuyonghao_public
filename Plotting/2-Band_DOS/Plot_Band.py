# coding: utf-8
#！/usr/bin/python3

'''
    plot bandstructure with vaspkit format(vaspkit.1.12)
    author: yonghao_zhu@163.com
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import linecache

################################################################################################################
'''
    INPUT NOTE:
    file_P = [];highkpt = [];color_bar_label = []
    file_P[0]:red  file_P[1]:blue(black)
    file_P[0]=color_bar_label[1]
    Subscript: '${MoS_2}$','${WS_2}$'
    Upscript: "DOS ($eV^{-1}$)"
'''
Plot_band_fun                    = {}
Plot_band_fun['TorF']            = True                    # band plot? True or False
Plot_band_fun['PlotMod']         = 'S'                     # 'C'->colorbar 'B'->bubble 'S'->line
#PlotMod=B 
Plot_band_fun['bubble_factor']   = 200                     # band_p0*bubble_factor for bubble mode
#PlotMod=C
Plot_band_fun['size_sca']        = 20                      # int
#General Settings
Plot_band_fun['show_']           = True                    # plot? True or False
Plot_band_fun['save_']           = False                   # save? True or False
Plot_band_fun['lw_']              = 2                      # line width
Plot_band_fun['y_range']         = [-2, 2, 1]          # [start, end, step]
Plot_band_fun['color_bar_label'] = ['${MoS_2}$','${WS_2}$']# NOTE: [blue,red]
Plot_band_fun['file_P']          = ['WS2.dat','MoS2.dat']  # NOTE: [red,blue]
Plot_band_fun['highkpt']         = ['M','K','Γ']       # same with KLABELS file, vaspkit format
Plot_band_fun['file_B']          = 'BAND.dat'              # vaspkit format
Plot_band_fun['file_K']          = 'KLABELS'               # vaspkit format
################################################################################################################
Plot_dos_fun                     = {}
Plot_dos_fun['TorF']             = False                   # dos plot? True or False
Plot_dos_fun['show_']            = True                    # plot? True or False
Plot_dos_fun['save_']            = False                   # save? True or False
Plot_dos_fun['Tdos_show']        = False                   # show Tdos? True or False
Plot_dos_fun['y_range']          = [0, 8, 2]               # [start, end, step]                      
Plot_dos_fun['x_range']          = [-3, 3, 1]              # [start, end, step]
Plot_dos_fun['lw_']              = 4                       # line width
Plot_dos_fun['file_T']           = 'TDOS.dat'              # TDOS file
Plot_dos_fun['file_PD']          = ['up.dat','dw.dat']     # PDOS file
Plot_dos_fun['legend_']          = ['${WS_2}$','${MoS_2}$']# legend name
Plot_dos_fun['color_PD']         = ['Red','Blue']          #PDOS line color
################################################################################################################

def readBand(file_B = 'BAND.dat', file_P = [], file_K = 'KLABELS'):
    try:
        print('---------------------------')
        line_2 = linecache.getline(file_B,2).strip('\n')
        print(line_2)
        nkpt,nband = int(line_2.split()[4]),int(line_2.split()[5])
        print('---------------------------')
        band_B = np.loadtxt(file_B)
        #test#print(band.shape)
        if len(file_P) >= 2:
            band_P = [np.loadtxt(i)[:,-1] for i in file_P]
        else:
            print('Two Ingredients!!!')
        p0 = band_P[0]/(band_P[0]+band_P[1])

        kl_file = open(file_K,'r');x_h = []
        for line in kl_file:
            lines = line.split()
            if len(lines) == 2:
                x_h.append(float(lines[1]))
        kl_file.close()
        print('high-kpt:','\n',x_h)
        print('---------------------------')
        return band_B, p0, x_h, band_P[0], band_P[1]
    except:
        print('ERROR in readBand!!!')

def readDos(file_T = 'TDOS.dat', file_PD = ['dw.dat','up.dat']):
    try:
        Tdos = np.loadtxt(file_T)
        Pdos = []
        if len(file_PD) >= 2:
            for i in file_PD:
                Pdos.append(np.loadtxt(i)[:,-1])
        return Tdos,Pdos
    except:
        print('ERROR in readDos!!!')

def PlotBand(y_range = [-3,3], show_ = False, size_sca = 20, save_ = False,
    file_B='BAND.dat',file_P=['WS2.dat','MoS2.dat'],file_K='KLABELS',
    highkpt=['Γ','M','K','Γ'],
    color_bar_label=['${MoS_2}$','${WS_2}$'],
    PlotMod='C',bubble_factor=60,lw_=2):

	##High Kpt of TMDs Materials(['Γ','M','K','Γ']), set highkpt for other materials.
	high_kpt = highkpt

	length = int((max(y_range[:2]) - min(y_range[:2])) / y_range[2])

	##Plot Band
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	if PlotMod == 'C':
		plt.figure(figsize=(7,5))
	else:
		plt.figure(figsize=(6,5))

	plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)

	#General Setting
	if PlotMod == 'C' or PlotMod == 'B':
		##read data
		band_B, band_P0, x_h, band_p0, band_p1 = readBand(file_B=file_B,file_P=file_P,file_K=file_K)
        #plot data, P0
		x = band_B[:,0];y = band_B[:,1];c = band_P0
		#plot band_line
		band_line = plt.plot(x,y,c='lightgrey',zorder=0,lw=lw_)
		#x and y range
		plt.xlim([min(x_h),max(x_h)]);plt.ylim(y_range[:2])
		#y label
		plt.ylabel("Energy (eV)",size=20)
		#x and y ticks(step=1)
		plt.xticks(x_h,high_kpt,size=20)

		plt.yticks([min(y_range[:2]) + i*y_range[2] for i in range(length+1)],
		    [min(y_range[:2]) + i*y_range[2] for i in range(length+1)],size=20)

		#high kpoint and fermi line
		x_h_ = x_h[1:-1]
		for i in x_h_:
		    plt.axvline(x=i,c="black",lw=1.5)
		plt.axhline(y=0,ls='--',c="black",lw=1,zorder=1)

	##color bar type 
	if PlotMod == 'C':       
		#color map: blue-red
		cmap = mpl.cm.jet
		#plot scatter with cmap and show color bar
		band_scatter = plt.scatter(x,y,c=c,s=size_sca,cmap=cmap)
		bar_color = plt.colorbar()
		#color bar set,
		bar_color.set_ticks([min(c),max(c)])
		bar_color.set_ticklabels(color_bar_label)
		bar_color.ax.tick_params(labelsize=14)
		#save .svg 
		if save_:
		    plt.savefig('band_color.png', dpi=600, format='png')

	##bubble type
	if PlotMod == 'B':
		#set legend
		
		##band bubble
		plt.scatter(x,y,s=bubble_factor*band_p0,color='r',alpha=0.55,linewidth=0, label=color_bar_label[1])
		plt.scatter(x,y,s=bubble_factor*band_p1,color='b',alpha=0.55,linewidth=0, label=color_bar_label[0])
		plt.legend(fontsize='large')
		#save .svg 
		if save_:
		    plt.savefig('band_bubble.png', dpi=600, format='png')

	##simple band type
	if PlotMod == 'S':
		band_B = np.loadtxt(file_B)
		kl_file = open(file_K,'r');x_h = []
		for line in kl_file:
		    lines = line.split()
		    if len(lines) == 2:
		        x_h.append(float(lines[1]))
		kl_file.close()
		#plot band_line
		x = band_B[:,0];y = band_B[:,1]
		plt.plot(x,y,c='black',zorder=0,lw=lw_)
		#x and y range
		plt.xlim([min(x_h),max(x_h)]);plt.ylim(y_range[:2])
		#y label
		plt.ylabel("Energy (eV)",size=20)
		#x and y ticks(step=1)
		plt.xticks(x_h,high_kpt,size=15)

		plt.yticks([min(y_range[:2]) + i*y_range[2] for i in range(length+1)],
		    [min(y_range[:2]) + i*y_range[2] for i in range(length+1)],size=20)

		#high kpoint and fermi line
		x_h_ = x_h[1:-1]
		for i in x_h_:
		    plt.axvline(x=i,c="black",lw=1.5)
		plt.axhline(y=0,ls='--',c="black",lw=1,zorder=1)
		##Set the thickness of the coordinate axis
		ax=plt.gca()
		ax.spines['bottom'].set_linewidth(1)
		ax.spines['left'].set_linewidth(1)
		ax.spines['right'].set_linewidth(1)
		ax.spines['top'].set_linewidth(1)
		#save .svg
		if save_:
		    plt.savefig('band_simple_band.png', dpi=600, format='png')

	if show_:
	    plt.show()

def PlotDos(y_range = [0,3], x_range = [0,8], 
    file_T = 'TDOS.dat', file_PD = ['dw.dat','up.dat'],
    lw_ = 2, legend_ = [],Tdos_show = False, color_PD = [], show_ = False, save_ = False):

    Tdos,Pdos = readDos(file_T=file_T,file_PD=file_PD)
    #test#print(Tdos.shape,len(Pdos))
    ##plot dos
    plt.rc('font',family='Times New Roman')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    plt.figure(figsize=(6.5,5))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.15)

    #plot Tdos
    x = Tdos[:,0];y_Tdos = Tdos[:,1]
    if Tdos_show == 'T':
        Tdos_line = plt.plot(x,y_Tdos,c='black',zorder=0,lw=lw_)
    #x and y range
    plt.xlim(x_range[:2]);plt.ylim(y_range[:2])
    #x and y ticks
    x_ = [min(x_range[:2])+i*x_range[2] for i in range(int((max(x_range[:2])-min(x_range[:2]))/x_range[2])+1)]
    y_ = [min(y_range[:2])+i*y_range[2] for i in range(int((max(y_range[:2])-min(y_range[:2]))/y_range[2])+1)]

    plt.xticks(x_,x_,size=15)
    plt.yticks(y_,y_,size=15)
    #x and y label
    plt.ylabel("DOS ($eV^{-1}$)",size=20)
    plt.xlabel("Energy (eV)",size=20)
    ##Plot PDOS
    for i in range(len(file_PD)):
        line_p = plt.plot(x,Pdos[i],c=color_PD[i],lw=lw_)
    #set legend
    plt.legend(legend_,loc=1,fontsize='large')#handles=[line_p],
    #save .svg 
    if save_:
        plt.savefig('dos.png', dpi=600, format='png')
    if show_:
        plt.show()

################################################################################################################
def main():
	##Plot_band
	if Plot_band_fun['TorF']:
	    PlotBand(PlotMod=Plot_band_fun['PlotMod'], y_range=Plot_band_fun['y_range'],
	        show_=Plot_band_fun['show_'], size_sca=Plot_band_fun['size_sca'], save_=Plot_band_fun['save_'],
	        file_P=Plot_band_fun['file_P'], file_B=Plot_band_fun['file_B'], file_K=Plot_band_fun['file_K'],
	        highkpt=Plot_band_fun['highkpt'], color_bar_label=Plot_band_fun['color_bar_label'],
	        bubble_factor=Plot_band_fun['bubble_factor'],lw_=Plot_band_fun['lw_'])

	##Plot_dos
	if Plot_dos_fun['TorF']:
	    PlotDos(y_range=Plot_dos_fun['y_range'], x_range=Plot_dos_fun['x_range'],
	        file_T=Plot_dos_fun['file_T'], file_PD=Plot_dos_fun['file_PD'],
	        lw_=Plot_dos_fun['lw_'], legend_=Plot_dos_fun['legend_'], color_PD=Plot_dos_fun['color_PD'],
	        save_=Plot_dos_fun['save_'], show_=Plot_dos_fun['show_'], Tdos_show=Plot_dos_fun['Tdos_show'])

if __name__ == '__main__':
	main()
