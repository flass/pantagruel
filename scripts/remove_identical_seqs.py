#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""remove duplicate sequences from a FASTA file based on a table of clusters of identical sequences"""
import sys

nfastain = sys.argv[1]
nfidentseq = sys.argv[2]
nfastaout = sys.argv[3]

lredundant = []
sfam = set()

if len(sys.argv)>4: 
	nfprevdbidentseq = sys.argv[4]
	print "will select representative sequences based on previous clustering"
	dprevfamreprseq = {}
	with open(nfprevdbidentseq, 'r') as fprevdbidentseq:
		for line in fprevdbidentseq:
			fam, prot = line.rstrip('\n').split('\t')
			if fam not in dprevfamreprseq:
				dprevfamreprseq[fam] = prot
	
	with open(nfidentseq, 'r') as fidentseq:
		for line in fidentseq:
			fam, prot = line.rstrip('\n').split('\t')
			if (fam in dprevfamreprseq):
				if prot != dprevfamreprseq[fam]:
					lredundant.append(prot)
			else:
				if fam not in sfam:
					sfam.add(fam)
				else:
					lredundant.append(prot)
else:
	with open(nfidentseq, 'r') as fidentseq:
		for line in fidentseq:
			fam, prot = line.rstrip('\n').split('\t')
			if fam not in sfam:
				sfam.add(fam)
			else:
				lredundant.append(prot)
			#~ prots = line.rstrip('\n').split('\t')
			#~ lredundant += prots[1:]

print "listed %d redundant sequences in dataset"%len(lredundant)
sredundant = set(lredundant)
print "generated hash index"
fastain = open(nfastain, 'r')
fastaout = open(nfastaout, 'w')
print "parsing redundant sequence fasta"
w = True
nnrseq = 0
for line in fastain:
	if line.startswith('>'):
		prot = line.strip('>\n').split(' ', 1)[0]
		if prot not in sredundant:
			w = True
			nnrseq += 1
		else:
			w = False
	if w:
		fastaout.write(line)

print "filtered %d non-redundant sequences"%nnrseq


fastain.close()
fastaout.close()
