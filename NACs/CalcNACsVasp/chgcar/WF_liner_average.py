import numpy as np
import linecache

# inputs
wf_or_charge = 'CHGCAR'
direction = 3 # not change, 3->z is ok
save_ = True
file_name = 'Liner' # not change

# number of atoms
line_7 = linecache.getline(wf_or_charge, 7)
atoms = sum([int(i) for i in line_7.split()])
print('atoms =', atoms)

# z, direction = 3 is ok 
line_5 = linecache.getline(wf_or_charge, 5)
z = float(line_5.split()[-1])

# grid 
line_grid = linecache.getline(wf_or_charge, 10+atoms)
grid = [int(i) for i in line_grid.split()]
print('grid =', grid)

# read grid data 
line_data = len(linecache.getline(wf_or_charge, 11+atoms).split())
print('line_data =', line_data)

tmp_ = ( grid[0] * grid[1] * grid[2] ) % line_data

if tmp_ == 0:
	number_of_lines = int( ( grid[0] * grid[1] * grid[2] ) / line_data )
else:
	number_of_lines = int( ( grid[0] * grid[1] * grid[2] ) / line_data ) + 1

print('total_lines =', number_of_lines)

# read
data_p = []

for num in range(number_of_lines):
	line = linecache.getline(wf_or_charge, 11+atoms+num).strip('\n').split()
	for ii in line:
		data_p.append(float(ii))

# check data	
if grid[0] * grid[1] * grid[2] != len(data_p):
	print('Error! grid[0] * grid[1] * grid[2] != len(data_p)')
	exit()
else:
	print('grid[0] * grid[1] * grid[2] = len(data_p)')

# split
data = np.array(data_p).reshape([grid[2], grid[1], grid[0]])

if save_:
	data_save = []
	name = file_name+'-'+wf_or_charge.replace('.vasp','.dat')
	with open(name,'w+') as f:
		for i in range(grid[direction-1]):
			if direction == 3:
				x = np.sum(data[i,:,:].reshape([1,-1]))
			elif direction == 2:
				x = np.sum(data[:,i,:].reshape([1,-1]))
			else:
				x = np.sum(data[:,:,i].reshape([1,-1]))

			data_save.append(x)

			f.writelines('  '+str(i*z/grid[direction-1])+'  '+str(x)+'\n')

	


