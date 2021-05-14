#! /usr/bin/python3
# --* coding: utf-8 *--

import numpy as np

high_kpoint = [[0.00, 0.00, 0.00],
               [0.00, 0.50, 0.00],
               [0.333333, 0.333333, 0.00],
			   [0.00, 0.00, 0.00]]

path_insert = 31

###################################
path = []
for k in range(len(high_kpoint)-1):
	path.append([high_kpoint[k], high_kpoint[k+1]])

kpoint = np.zeros([(len(high_kpoint)-1)*(2+path_insert),4])

k_list = []
for i in path:
	point_1 = np.array(i[0])
	point_2 = np.array(i[1])
	delt = (point_2 - point_1)/(path_insert+1)
	for i in range(path_insert+2):
		point = point_1 + delt*i
		k_list.append(point)

for i in range(len(k_list)):
	kpoint[i,:3] = k_list[i]

np.savetxt('add.kpt', kpoint, fmt='%.06f')

