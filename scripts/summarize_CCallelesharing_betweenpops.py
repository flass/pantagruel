#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys, os

nflnfpopref = sys.argv[1]
nfpops = sys.argv[2]
nfout = sys.argv[3]

with open(nflnfpopref, 'r') as flnfpopref:
	lnfpopref = [line.rstrip('\n') for line in flnfpopref]

lpops = []
with open(nfpops, 'r') as fpops:
	for line in fpops:
		if line.startswith('#'): continue
		lpops.append(line.split('\t')[0])

fout = open(nfout, 'w')
fout.write('\t'.join(['']+lpops)+'\n')

for nfpopref in lnfpopref:
	fam = os.path.basename(nfpopref).split('-')[0]
	with open(nfpopref, 'r') as fpopref:
		for line in fpopref:
			colclade, spopfreqs = line.rstrip('\n').split('\t', 1)
			famcc = fam+'-'+colclade
			dpopfreqs = {p:float(f) for p, f in [pf.split(':') for pf in spopfreqs.split('\t')]}
			fout.write('\t'.join([famcc]+[str(dpopfreqs.get(pop, '0.0')) for pop in lpops])+'\n')

fout.close()
		
		
