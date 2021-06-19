# coding: utf-8
#！/usr/bin/python3

'''
    plot bandstructure with vaspkit format(vaspkit.1.12)
    author: yonghao_zhu@163.com
'''

import os

#######################################################
SplitBand_fun          = {}
SplitBand_fun['file']  = 'BAND-spin.dat'
#######################################################

def SplitBand(file = 'BAND-spin.dat'):
	'''
	out: band_up.dat and band_dw.dat
	BAND.dat format:
	#K-Path(1/A)         Spin-Up(eV)   Spin-Down(eV)
	# NKPTS & NBANDS:  93  32
	# Band-Index    1
	# Band-Index    2
	'''

	spin_up,spin_dw = [],[]

	with open(file, 'r') as band:
		for line in band:
			lines = line.split()
			if '#' in lines or 'Spin-Up(eV)' in lines:
				spin_dw.append(line)
				spin_up.append(line)
			elif len(lines) == 3:
				up = '   ' + lines[0] + '    '+lines[1] + '\n'
				dw = '   ' + lines[0] + '    '+lines[2] + '\n'
				spin_up.append(up)
				spin_dw.append(dw)

	with open('band_up.dat', 'w+') as f:
		for i in spin_up:
			f.writelines(i)

	with open('band_dw.dat', 'w+') as f:
		for i in spin_dw:
			f.writelines(i)

def main():
	if os.path.isfile(SplitBand_fun['file']):
		SplitBand(file = SplitBand_fun['file'])
	else:
		print('No File!')

if __name__ == '__main__':
	main()