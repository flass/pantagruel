#!/usr/bin/python
import sys, glob
from os import path

target = sys.argv[1]

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

lnf = glob.glob(target)
nfout = '_'.join(target.rsplit('/', 1)).replace('*', '_').rstrip('_')+'_line_counts'
fout = open(nfout, 'w')

for nf in lnf:
	with open(nf, 'r') as f:
		fout.write(['\t'.join((path.basename(nf), str(file_len(nf))))+'\n')

fout.close()

		
