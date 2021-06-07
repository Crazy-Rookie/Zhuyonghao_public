#! /usr/bin/python3
# -*- conding=UTF-8 -*-
#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /  Zhu  \
#    //  Yong \\
#   //|  Hao  |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__)))  

'''
	1. plot sd and dephasing
	2. yonghao_zhu@163.com
'''

#################################################
import os                                       # 
import numpy as np                              #
from scipy import optimize                      #
from scipy.constants import physical_constants  #
import matplotlib as mpl                        #
import matplotlib.pyplot as plt                 #
#################################################

#################################################
paras                 = {}
paras['energy_file']  = 'energ1'
paras['ene_ini_fin']  = [1,2] #start 1, including  
#paras['time_fin'] - paras['time_ini'] > 800 
paras['time_ini']     = 1     #start 1, including
paras['time_fin']     = 1500  #start 1, including
paras['mode']         = 'sd' # un_acf, dephasing, sd
show_                 = True  # True or False
save_                 = False # True or False
#################################################

def AutoCorrectionFunction(paras, save_ = False, show_ = False):
	'''
	delta_Eij = Eij - <Eij>
	Cij = <delta_Eij(t)delta_Eij(t-t')>
	'''

	ene_p = np.loadtxt(paras['energy_file'])

	if len(paras['ene_ini_fin']) == 1:
		gap = ene_p[paras['time_ini']-1: paras['time_fin']]

	elif len(paras['ene_ini_fin']) == 2:
		gap_1 = ene_p[paras['time_ini']-1: paras['time_fin'], paras['ene_ini_fin'][0]-1]
		gap_2 = ene_p[paras['time_ini']-1: paras['time_fin'], paras['ene_ini_fin'][1]-1]

		gap = (gap_2 - gap_1).T

	else:
		print('Check paras[\'ene_ini_fin\']')
		exit()
	
	delta_gap = gap - np.mean(gap)

	times = paras['time_fin'] - paras['time_ini'] + 1

	un_acf = np.correlate(delta_gap, delta_gap, "full")[-times:] / (times) 
	
	acf = un_acf / un_acf[0]

	if show_:

		x1 = [i for i in range(0, 800)]
		y2 = [0 for i in range(800)]

		plt.figure(figsize=(7,5))
		plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)
		
		plt.xlabel('Time (fs)', fontsize=16)
		plt.ylabel('C$_u$$_n$(t) (eV$^2$)', fontsize=16)
		
		plt.xlim(0, 800)

		plt.plot(x1, un_acf[:800], 'r-')
		plt.plot(x1, y2, 'r--')
		
		plt.show()

	if save_:

		np.savetxt(paras['mode']+'.txt', un_acf, fmt='%.6f')
	
	return un_acf, acf, gap

#################################################

def gauss_function(x, a):
    return np.exp(-0.5 * (x / a) ** 2)

def Depahsing(paras, save_ = False, show_ = False):

	hbar = 1e15 * physical_constants['Planck constant over 2 pi in eV s'][0]

	print('hbar=', hbar)
	print('linewidth=hbar/dephasing-time', )

	un_acf = AutoCorrectionFunction(paras)[0]

	integral2 = np.array([np.sum(un_acf[0:i]) for i in range(np.size(un_acf))]) / hbar

	integral1 =np.array([np.sum(integral2[0:i]) for i in range(np.size(un_acf))]) / hbar
	
	deph = np.exp(-integral1)

	if show_:

		plt.figure(figsize=(7,5))
		plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

		plt.xlabel('Time (fs)', fontsize=16)
		plt.ylabel('Dephasing', fontsize=16)
		
		time = 200

		plt.xlim(0, time)  
		plt.ylim(0, 1)  

		popt, pcov = optimize.curve_fit(
		    gauss_function, [i for i in range(0, time)], deph[1:time + 1])

		plt.title('dephasing: %.4f fs' % popt[0] + ',  linewidth: %.6f' %(hbar / popt[0]))

		plt.plot([i for i in range(0, time + 2)], deph[1:time + 1 + 2], 'r-')
		plt.show()

	if save_:
		np.savetxt(paras['mode']+'.txt', deph, fmt='%.6f')

#################################################

def SD(paras, save_ = False, show_ = False):

	un_acf, acf, gap = AutoCorrectionFunction(paras, save_ = False, show_ = False)

	sd_y = abs(1 / np.sqrt(2 * np.pi) * np.fft.fft(un_acf, 100000)) ** 2
	
	sd_x = np.fft.fftfreq(len(sd_y))
	
	sd_x = sd_x * 33356.40952  # fs --> cm-1

	if show_:

		plt.figure(figsize=(7,5))
		plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

		plt.xlabel('Frequency (cm$^-$$^1$)', fontsize=16)
		plt.ylabel('Spectral Density', fontsize=16)
		
		plt.xlim(0, 800)  # x显示区间
		plt.ylim(0, max(sd_y))  # y显示区间
		
		plt.plot(sd_x, sd_y, 'r-')
		
		plt.show()

	if save_:
		sd = np.stack([sd_x, sd_y]).T[:4000,:]
		np.savetxt(paras['mode']+'.txt', sd, fmt='%.6f')
	

#################################################
def main():

	if not os.path.isfile(paras['energy_file']):
		print('No paras[\'energy_file\']')
		exit()


	if paras['mode'] == 'un_acf':
		AutoCorrectionFunction(paras, 
							   save_ = save_,
							   show_ = show_)

	if paras['mode'] == 'dephasing':
		Depahsing(paras,
			      save_ = save_,
			      show_ = show_)

	if paras['mode'] == 'sd':
		SD(paras,
		   save_ = save_,
		   show_ = show_)


if __name__ == '__main__':
	main()