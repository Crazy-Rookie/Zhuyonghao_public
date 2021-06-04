# coding: utf-8
#！/usr/bin/python3
'''
	1. author: yonghao_zhu@163.com
	2. plot time-xxx, such as time-NAC
	3. For NAC-time plot, only Line style until 2020.10.06	
'''
##########################################
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
##########################################
##########################################
PlotNACTime_fun                   = {}
PlotNACTime_fun['TorF']           = True         # True or False
PlotNACTime_fun['path_']          = './'
PlotNACTime_fun['file_']          = ['NAC_Time-NormalRun.dat','NAC_Time-PhaseCorrRun.dat']
PlotNACTime_fun['Scaling_factor'] = [1.00, 1.00] #0.33 or 1
PlotNACTime_fun['legend_']        = ['Normal','Phase Corr.']
PlotNACTime_fun['color']          = ['black','blue']
PlotNACTime_fun['PlotMod']        = 'Line'       #unitl now only Line style
#PlotMod=Line
PlotNACTime_fun['Time_scale']     = [3000,5000]  #Include boundary(x range)
PlotNACTime_fun['Absolute_value'] = 'F'    #T-> absolute value
PlotNACTime_fun['NAC_scale']      = [[-4.5,4.5],[-1.5,1.5]]#meV Include boundary(y range)
PlotNACTime_fun['show_']          = True   # True or False
PlotNACTime_fun['save_']          = False  # True or False
##########################################
##########################################
PlotBandEneTime_fun               = {}
PlotBandEneTime_fun['TorF']       = True   # True or False
PlotBandEneTime_fun['file_']      = 'energy_new'
PlotBandEneTime_fun['Time_scale'] = [1,1500]
PlotBandEneTime_fun['Ene_scale']  = [-1.8,-0.8]
PlotBandEneTime_fun['show_']      = True   # True or False
PlotBandEneTime_fun['save_']      = False  # True or False
##########################################
##########################################
PlotTotEneTime_fun               = {}
PlotTotEneTime_fun['TorF']       = True   # True or False
PlotTotEneTime_fun['file_']      = 'TimeTemEn.dat'
PlotTotEneTime_fun['Time_scale'] = [1,3000]
PlotTotEneTime_fun['show_']      = True   # True or False
PlotTotEneTime_fun['save_']      = False  # True or False
PlotTotEneTime_fun['TemporEne']  = 'temp'  #temp->temperature-time or ene->energy-time
##########################################
##########################################
PlotIPRTime_fun                  = {}
PlotIPRTime_fun['TorF']          = True  # True or False
PlotIPRTime_fun['path_']         = './ipr/'
PlotIPRTime_fun['file_perfix']   = 'ipr'
#time in [0,9999]
PlotIPRTime_fun['ini_time'] = 2 #in
PlotIPRTime_fun['fin_time'] = 2000 #in
PlotIPRTime_fun['ini_band'] = 20 #in(for example: 70-50=20)
PlotIPRTime_fun['fin_band'] = 21 #in(for exampel: 71-50=21)
PlotIPRTime_fun['legend_']  = ['VBM', 'CBM']
PlotIPRTime_fun['show_']    = True  # True or False
PlotIPRTime_fun['y_scale']  = [0,0.001]
PlotIPRTime_fun['save_']    = False  # True or False
##########################################
##########################################
PlotOthers_fun              = {}
PlotOthers_fun['TorF'] = 'F'
##########################################
def PlotNACTime(path_='./', file_=['NAC-Time.dat'], PlotMod='Line', Time_scale=[2,2000],
	Absolute_value=False, NAC_scale=[-3,3], show_=True, save_=False,
	x_range=[1,2], y_range=[-1,1], Scaling_factor=[0.33], legend_=['Phase Corr'],
	color=[]):

	file = [path_+i for i in file_]

	if PlotMod == 'Line':

		plt.rc('font',family='Times New Roman')
		mpl.rcParams['xtick.direction'] = 'in'
		mpl.rcParams['ytick.direction'] = 'in'

		if len(file) == 1:#one nac-time plot

			#read file
			nac_data = np.loadtxt(file[0])
			#test#print(nac_data.shape)
			time_x = nac_data[:,0]
			if Absolute_value:
				nac_y = abs(nac_data[:,1])*Scaling_factor[0]*1000#ev->meV
			else:
				nac_y = nac_data[:,1]*Scaling_factor[0]*1000#ev->meV

			plt.figure(figsize=(7,5))
			plt.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.12)

			##plot nac-time
			plt.plot(time_x,nac_y,c='black',zorder=0,lw=2)
			#set legend
			plt.legend(legend_,loc=1,fontsize='large')
			#xlim and ylim and xticks and yticks
			plt.xticks(size=15); plt.yticks(size=15)
			plt.xlim(Time_scale); plt.ylim(NAC_scale[0])
			#x and y label
			plt.ylabel("NAC (meV)",size=20); plt.xlabel("Time (fs)",size=20)

		else:#more than one nac-time plot

			num_nac = len(file)
			mpl.rcParams['axes.unicode_minus'] = False 
			fig,ax_big = plt.subplots(sharex=True)
			axis = plt.gca()
			plt.subplots_adjust(wspace=0, hspace=0)

			axis.spines['right'].set_color('none')
			axis.spines['top'].set_color('none')
			axis.spines['bottom'].set_color('none')
			axis.spines['left'].set_color('none')

			axis.spines['bottom'].set_position(('outward', 30))
			axis.spines['left'].set_position(('outward', 30))

			ax_big.set_xticks([]) 
			ax_big.set_yticks([])

			ax_big.set_ylabel("NAC (meV)",size=20)
			plt.xlabel("Time (fs)",size=20)

			plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
			
			for i in range(1, 1 + num_nac):
				nac_data = np.loadtxt(file[i-1]); time_x = nac_data[:,0]			
				if Absolute_value:
					nac_y = abs(nac_data[:,1])*Scaling_factor[i-1]*1000#ev->meV
				else:
					nac_y = nac_data[:,1]*Scaling_factor[i-1]*1000#ev->meV
				ax = fig.add_subplot(num_nac, 1, i)
				ax.plot(time_x,nac_y,c=color[i-1])
				ax.set_ylim(NAC_scale[i-1])
				if i == num_nac:
					ax.set_xlim(Time_scale)
				else:
					ax.set_xlim(Time_scale)
					ax.set_xticks([]) 
				#set legend
				ax.legend((legend_[i-1],),loc=1,fontsize='large')

		#save .svg
		if save_:
			plt.savefig('NAC-Time.svg',format='svg')
		if show_:
			plt.show()

	else:
		print('Only Line type Now!!!')
	print(1)

if PlotNACTime_fun['TorF']:
	if os.path.isfile(PlotNACTime_fun['file_'][0]):
		PlotNACTime(path_=PlotNACTime_fun['path_'], file_=PlotNACTime_fun['file_'], PlotMod=PlotNACTime_fun['PlotMod'],
			Time_scale=PlotNACTime_fun['Time_scale'], Absolute_value=PlotNACTime_fun['Absolute_value'],
			NAC_scale=PlotNACTime_fun['NAC_scale'], show_=PlotNACTime_fun['show_'], save_=PlotNACTime_fun['save_'],
			Scaling_factor=PlotNACTime_fun['Scaling_factor'], legend_=PlotNACTime_fun['legend_'], color=PlotNACTime_fun['color'])
	else:
		print('No %s!' % PlotNACTime_fun['file_'][0])

#********************************************************************************************
#********************************************************************************************
def PlotBandEneTime(file_='x', Time_scale=[1,2], Ene_scale=[1,2], show_=True, save_=False):
	energy_data = np.loadtxt(file_)
	#test#print(energy_data.shape)
	time = energy_data.shape[0]
	states = energy_data.shape[1]
	time_x = np.array([i+1 for i in range(time)])
	#plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'
	plt.plot(time_x,energy_data,lw=2)

	plt.figure(figsize=(7,5))
	plt.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.12)

	#xlim and ylim and xticks and yticks
	plt.xticks(size=15); plt.yticks(size=15)
	plt.xlim(Time_scale); plt.ylim(Ene_scale[0])
	#x and y label
	plt.ylabel("Energy (eV)",size=20); plt.xlabel("Time (fs)",size=20)

	#save .svg
	if save_:
		plt.savefig('NAC-Time.svg',format='svg')
	if show_:
		plt.show()

if PlotBandEneTime_fun['TorF']:
	if os.path.isfile(PlotBandEneTime_fun['file_']):
		PlotBandEneTime(file_=PlotBandEneTime_fun['file_'],
			Time_scale=PlotBandEneTime_fun['Time_scale'],
			Ene_scale=PlotBandEneTime_fun['Ene_scale'],
			show_=PlotBandEneTime_fun['show_'], save_=PlotBandEneTime_fun['save_'])
	else:
		print('No %s!' % PlotBandEneTime_fun['file_'])
#********************************************************************************************
#********************************************************************************************

def PlotTotEneTime(file_='TimeTemEn.dat', Time_scale=[1,2], show_=True, save_=False, TemporEne='temp'):
	data = np.loadtxt(file_,comments=['Time'])
	#test#print(data.shape)
	time_x = data[:,0]
	if TemporEne == 'temp':#temperature-time plot
		data_y = data[:,2]
	if TemporEne == 'ene':#energy-time plot
		data_y = data[:,1]
	#plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	plt.figure(figsize=(7,5))
	plt.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.12)

	plt.plot(time_x,data_y,c='black',lw=2)
	#xlim and ylim and xticks and yticks
	plt.xticks(size=15); plt.yticks(size=15)
	#x and y label
	if TemporEne == 'temp':#temperature-time plot
		plt.ylabel("Temperature (K)",size=20); plt.xlabel("Time (fs)",size=20)
		plt.xlim(Time_scale); plt.ylim([0,600])
	if TemporEne == 'ene':#energy-time plot
		plt.ylabel("Energy (eV)",size=20); plt.xlabel("Time (fs)",size=20)
		plt.xlim(Time_scale); plt.ylim([float(data_y.min())-0.5,float(data_y.max())+0.5])

	#save .svg
	if save_:
		plt.savefig('NAC-Time.svg',format='svg')
	if show_:
		plt.show()

if PlotTotEneTime_fun['TorF']:
	if os.path.isfile(PlotTotEneTime_fun['file_']):
		PlotTotEneTime(file_=PlotTotEneTime_fun['file_'], Time_scale=PlotTotEneTime_fun['Time_scale'],
			show_=PlotTotEneTime_fun['show_'], save_=PlotTotEneTime_fun['save_'],
			TemporEne=PlotTotEneTime_fun['TemporEne'])
	else:
		print('No %s!' %PlotTotEneTime_fun['file_'])
#********************************************************************************************
#********************************************************************************************

def PlotIPRTime(path_='./ipr/', ini_time=2, fin_time=3, file_perfix='ipr',
				show_=True, save_=False, ini_band=20, fin_band=21, y_scale=[0,0.0002],
				legend_=['1','2']):
	#read files
	ipr_data = np.zeros((fin_time-ini_time+1, fin_band-ini_band+1))
	#test#
	print(ipr_data.shape)
	for time in range(ini_time,fin_time+1):
		ipr_name = path_+file_perfix+'%.04i' % time
		'''
		band_ini
		band_ini+1
		band_ini+2
		...
		'''
		ipr_ = np.loadtxt(ipr_name)[ini_band:fin_band+1].reshape(fin_band-ini_band+1,-1).T
		ipr_data[time-ini_time,:] = ipr_
	#test#print(ipr_data[0])
	
	plt.figure(figsize=(7,5))
	plt.subplots_adjust(left=0.18, right=0.9, top=0.9, bottom=0.12)

	#plot
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'
	time_x = [i for i in range(ini_time, fin_time+1)]
	plt.plot(time_x,ipr_data,lw=2)
	#xlim and ylim and xticks and yticks
	plt.xticks(size=15); plt.yticks(size=15)
	plt.xlim([ini_time, fin_time]); plt.ylim(y_scale)
	#legend
	plt.legend(legend_,loc=1,fontsize='large')
	#x and y label
	plt.ylabel("IPR",size=20); plt.xlabel("Time (fs)",size=20)
	plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
	#save .svg
	if save_:
		plt.savefig('IPR-Time.png',format='png',dpi=300)
	if show_:
		plt.show()

if PlotIPRTime_fun['TorF']:
	if os.path.isdir(PlotIPRTime_fun['path_']):
		PlotIPRTime(path_ = PlotIPRTime_fun['path_'],
					file_perfix = PlotIPRTime_fun['file_perfix'],
					ini_time = PlotIPRTime_fun['ini_time'],
					fin_time = PlotIPRTime_fun['fin_time'],
					show_ = PlotIPRTime_fun['show_'],
					save_ = PlotIPRTime_fun['save_'],
					ini_band = PlotIPRTime_fun['ini_band'],
					fin_band = PlotIPRTime_fun['fin_band'],
					y_scale = PlotIPRTime_fun['y_scale'],
					legend_ = PlotIPRTime_fun['legend_'])

	else:
		print('No %s!' % PlotIPRTime_fun['path_'])

#********************************************************************************************
#********************************************************************************************
def PlotOthers():
	print('Comimg Soon ...')	

if PlotOthers_fun['TorF']:
	PlotOthers()