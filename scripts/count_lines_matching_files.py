import sys, glob
from os import path

target = sys.argv[1]

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

lnf = glob.glob(target)

lnumlines = []

for nf in lnf:
	with open(nf, 'r') as f:
		lnumlines.append((path.basename(nf), file_len(nf)))

with open(target.replace('*', '_')+'line_counts', 'w') as fout:
	fout.write('\n'.join(['\t'.join(t) for t in lnumlines]))


		
