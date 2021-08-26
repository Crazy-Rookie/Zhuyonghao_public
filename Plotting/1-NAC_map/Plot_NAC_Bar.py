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

##################################
import numpy as np               #
import matplotlib as mpl         #
import matplotlib.pyplot as plt  #
import os                        #
from matplotlib import cm        #
##################################

##################################
params = {}
# files: two columns (time NAC), set start time to 1
params['file_name']  = ['NAC_Time_Corr_POT.dat', 'NAC_Time_HSE06.dat', 'NAC_Time_Normal_POT.dat']
# legend
params['legend']     = ['GGA-1/2', 'HSE06', 'PBE']
params['time_range'] = [1, 1400, 400] # [start, end, interval] including, start 1
# color bar range
params['colorbar_range'] = [0, 4, 5] # [minmum, maxmum, numbers] meV
params['show_']      = True  # True or False
params['save_']      = False

##################################
def PlotNACBar(params):
	'''
		plot nac (absolute value) to bar type, like Nano Lett. 2021, 21, 10, 4403–4409.
	'''

	# read data
	nac_data = params['nac_data']

	# normalization of NAC value
	norm = plt.Normalize(0, params['colorbar_range'][1])

	# time 
	time_x = np.array([i for i in range(params['time_range'][0], params['time_range'][1]+1)])
	#test#print(time_x)

	colorbar_range = params['colorbar_range']

	# y
	y = np.array([params['colorbar_range'][1] for i in range(params['time_range'][0], params['time_range'][1]+1)])

	# plot setting
	plt.rc('font',family='Times New Roman')
	mpl.rcParams['xtick.direction'] = 'in'
	mpl.rcParams['ytick.direction'] = 'in'

	# one plot
	if len(nac_data) == 1:

		plt.figure(figsize=(7,5))
		plt.subplots_adjust(left=0.1, right=0.9, top=0.8, bottom=0.2)
		
		# NAC: absolute value
		nac_y = nac_data[0]
		nac_y = norm(nac_y)
		#test#print(nac_y)

		# plot: color bar
		map_vir = cm.get_cmap(name='Reds')
		colors = map_vir(nac_y)
		plt.bar(time_x, y, color=colors, width=1)
		sm = cm.ScalarMappable(cmap=map_vir, norm=norm)
		
		cbar = plt.colorbar(sm, ticks=np.linspace(colorbar_range[0], colorbar_range[1], colorbar_range[2]))
		cbar.set_label('NAC (meV)', fontsize=20)
		cbar.ax.tick_params(labelsize=18)

		# xlim and ylim
		plt.xlim(params['time_range'][:2]); plt.ylim(params['colorbar_range'][:2])
		# xticks and yticks
		plt.yticks([])

		plt.xlabel("Time (fs)",size=25)

		length = int((max(params['time_range'][:2]) - min(params['time_range'][:2])) / params['time_range'][2])
		x_range = [params['time_range'][0]-1, params['time_range'][1]-1]
		plt.xticks([min(x_range) + i*params['time_range'][2] for i in range(length+1)],
			       [min(x_range) + i*params['time_range'][2] for i in range(length+1)], size=20)

	# more than one plot
	if len(nac_data) != 1:		
		
		num_nac = len(nac_data)

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

		plt.xlabel("Time (fs)",size=25)

		plt.subplots_adjust(left=0.1, right=0.84, top=0.93, bottom=0.17)

		for i in range(1, 1 + num_nac):		

			ax = fig.add_subplot(num_nac, 1, i)
			
			nac_y = nac_data[i-1]

			nac_y = norm(nac_y)

			map_vir = cm.get_cmap(name='Reds')
			colors = map_vir(nac_y)

			plt.bar(time_x, y, color=colors, width=1)

			ax.set_yticks([])
			ax.set_ylim(params['colorbar_range'][:2])
			ax.set_xlim(params['time_range'][:2])

			if i == num_nac:
				length = int((max(params['time_range'][:2]) - min(params['time_range'][:2])) / params['time_range'][2])
				x_range = [params['time_range'][0]-1, params['time_range'][1]-1]
				plt.xticks([min(x_range) + i*params['time_range'][2] for i in range(length+1)],
					       [min(x_range) + i*params['time_range'][2] for i in range(length+1)], size=20)
			else:
				ax.set_xticks([])

		map_vir = cm.get_cmap(name='Reds')
		sm = cm.ScalarMappable(cmap=map_vir, norm=norm)

		# color bar location
		cb_location = [0.86, 0.2, 0.03, 0.6] # [left, bottom, width, hight]
		cbar_ax = fig.add_axes(cb_location)
		
		cbar = plt.colorbar(sm, cax=cbar_ax, ticks=np.linspace(colorbar_range[0], colorbar_range[1], colorbar_range[2]))
		cbar.set_label('NAC (meV)', fontsize=20)
		cbar.ax.tick_params(labelsize=18)

	if params['save_']:
		plt.savefig("NACMapBar.png", format='png', dpi=600)
	if params['show_']:
		plt.show()		

##################################
def main(params):

	times_inp = params['time_range'][1] - params['time_range'][0] + 1

	params['nac_data'] = []

	# check and read files
	for name in params['file_name']:
		# check file name
		if name not in os.listdir('./'):
			print(name)
			print('No file: %s! Exiting...' % name)
			exit()
	
		# check file data
		else:
			nac = np.loadtxt(name)
			times_read = nac.shape[0]
			if times_read < times_inp:
				print('The input times (%s) do not match the %s (%s)! Exiting...' %(params['time_range'], name, times_read))
				exit()
			elif params['time_range'][1] > times_read:
				print('The maximum of input times exceed the %s (%s)! Exiting...' %(name, times_read))
				exit()
			elif params['time_range'][0] > 0 and params['time_range'][1] <= times_read and  times_read >= times_inp:
				print('Reading (from eV to meV) the %s...' % name)
				params['nac_data'].append(np.abs(nac[:, 1])*1000)
			else:
				print('ERROR! Exiting...')
				exit()


	# normalization of NAC value
	print('-----------------------------------------')
	print('Normalization of NAC value...')
	max_colorbar_range = params['colorbar_range'][1]
	max_nac_read = max([np.max(i) for i in params['nac_data']])
	print('max_colorbar_range=', max_colorbar_range, 'meV')
	print('max_nac_read=', max_nac_read, 'meV')

	if max_colorbar_range < max_nac_read:
		print('Ensure max_colorbar_range >= max_nac_read! Exiting...')
		exit()

	# plot
	print('-----------------------------------------')

	PlotNACBar(params)

if __name__ == "__main__":
	main(params)

