# !/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import nrrd
import bonndit.directions.csd_peaks as cs


def main():
	parser = argparse.ArgumentParser(
		description='This script evaluates a 8th order tensor fodf on a sphere and selects k maxima with minimum seperation'
					' angle and minimum difference to the minimum. Then each peak is refined by a gradient descent. ')
	parser.add_argument('infile',
						help='4D input file containing fODFs in masked 8-order tensor format (1+45 fODF coefficients,x,y,z)')

	parser.add_argument('outfile',
						help='5D output file with the approximation result (4,r,x,y,z)')
	parser.add_argument('-p', help='peaks')
	parser.add_argument('-sa', help='Seperation angle in degrees')
	parser.add_argument('-h', help='minimum height')
	args = parser.parse_args()
	# Load fODF input file
	fodfs, meta = nrrd.read(args.infile)

	if fodfs.shape[0] != 46:
		raise Exception("fodf has to be 4th order tensor.")
	if len(fodfs.shape) != 4:
		raise Exception("fodf have to be in 3d space. Hence, fodf has to be 4d.")

	data = fodfs.reshape((fodfs.shape[0], -1))
	output = cs.csd_peaks(data, int(args.p), float(args.sa), float(args.h))

	output = output.reshape((4, int(args.p)) + fodfs.shape[1:])
	newmeta = {k: meta[k] for k in ['space', 'space origin']}
	newmeta['kinds'] = ['list', 'list', 'space', 'space', 'space']
	newmeta['space directions'] = np.vstack(([np.nan, np.nan, np.nan], meta['space directions']))
	nrrd.write(args.outfile, np.float32(output), newmeta)


if __name__ == "__main__":
	main()
