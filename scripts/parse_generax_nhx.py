#!/usr/bin/env python2.7

import tree2

import sys, os, getopt
from parseGeneRaxRec import parseMultiNHXRec

#sys.setrecursionlimit(20000)

if __name__ == "__main__":

	nfinnhx = sys.argv[1]
	nfout = sys.argv[2]

	devtcount, lgt = parseMultiNHXRec(nfinnhx)

	with open(nfout, 'w') as fout:
		for evt, n in devtcount.iteritems():
			fout.write('\t'.join(evt+(str(n),))+'\n')
