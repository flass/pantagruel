#!/usr/bin/env python2.7

import tree2

import sys, os, getopt
from parseGeneRaxRec import parseMultiNHXRec

def usage():
	s =  "Usage: %s -G rec_gene_tree_file.nhx -S reference_species_tree_file -o output_file"
	s += "options:"
	s += "		-g	will 'guess' infer additional HGT event from the pattern of annotation of the reconciled NHX gene tree file where reporting of events may have been overlooked (experimental!)"
	return s

#sys.setrecursionlimit(20000)

if __name__ == "__main__":

	opts, args = getopt.getopt(sys.argv[1:], 'G:S:o:gh', ['help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	nfinnhx = dopt['-G']
	nfspetree = dopt['-S']
	nfout = dopt['-o']
	guessHGT = ('-g' in dopt)
	
	spetree = tree2.AnnotatedNode(file=nfspetree, namesAsNum=True)
	#~ spetree = tree2.AnnotatedNode(file=nfspetree)
	#~ print spetree
	devtcount, lgt = parseMultiNHXRec(nfinnhx, spetree, guessHGT=guessHGT)

	with open(nfout, 'w') as fout:
		for evt, n in devtcount.iteritems():
			fout.write('\t'.join(evt+(str(n),))+'\n')
