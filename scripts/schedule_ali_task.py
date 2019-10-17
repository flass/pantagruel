#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os

nffamtab = sys.argv[1]
dirfamfasta = sys.argv[2]
nfoutalitasklit = sys.argv[3]
nparalleltasks = int(sys.argv[4])
excludefams = sys.argv[5:]	# ORFan families
fastaext = '.fasta'

# 2 nr sequences can translte into many sequqnces in whole genome db, so must have sequqnces aligned
minseq2ali = 2
dfamsize = {}

# lists (nr) gene family sizes
with open(nffamtab, 'r') as ffamtab:
	for line in ffamtab:
		lsp = line.rstrip('\n').split('\t')
		dfamsize[lsp[0]] = dfamsize.get(lsp[0], 0) + 1

famsizesummary = open("%s.sizes"%nffamtab, 'w')
# order families by size
dsizefam = {}
for fam in dfamsize:
	s = dfamsize[fam]
	famsizesummary.write("%s\t%d\n"%(fam, s))
	if fam in excludefams or s<minseq2ali: continue
	dsizefam.setdefault(s, []).append(fam)

famsizesummary.close()
	
famsizes = dsizefam.keys()
famsizes.sort(reverse=True)

dtaskout = {}
if nparalleltasks > 1:
	for t in range(nparalleltasks):
		dtaskout[t] = open("%s.%d"%(nfoutalitasklit, t), 'w')
else:
	dtaskout[0] = open(nfoutalitasklit, 'w')

# distribute task among task lists for similar-load parallel execution
n = 0
for famsize in famsizes:
	for fam in dsizefam[famsize]:
		dtaskout[n].write("%s/%s%s\n"%(dirfamfasta, fam, fastaext))
		n = (n + 1)%nparalleltasks

for t in dtaskout:
	dtaskout[t].close()
