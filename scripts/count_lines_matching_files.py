#!/usr/bin/python
import sys, glob
from os import path

target = sys.argv[1]

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

nfout = ('_'.join([path.dirname(target), path.basename(target)]).replace('*', '_')+'_line_counts').replace('__', '_')
print nfout
fout = open(nfout, 'w')
lnf = glob.glob(target)

for nf in lnf:
	with open(nf, 'r') as f:
		fout.write('\t'.join((path.basename(nf), str(file_len(nf))))+'\n')

fout.close()

		
