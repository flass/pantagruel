#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, getopt, os
from ptg_utils import parseChain
import multiprocessing as mp


def usage():
	s =  'Basic usage:\n'
	s += ' python %s -G /path/to/list.of.gene.tree.sample.chain1.files -o /path/to/output.folder [options]\n'%sys.argv[0]
	s =  'Options:\n'
	s += '  --chain.ext\t\tdefine the file extension of the tree chain file (the #1 chain when multiple ones); default to \'mb.run1.t\'\n'
	s += '  --nb.chains\t\tthe number of available tree chain to interwine into the single result file; default 2\n'
	s += '  --threads\t\t\tnumber of parralel processes to run; default to 1 (sequential).\n'
	s += '  -v  --verbose\tverbose mode.\n'
	return s

opts, args = getopt.getopt(sys.argv[1:], 'G:o:hv', ['--nb.chains=', '--chain.ext=', 'verbose=', 'help', 'threads='])
dopt = dict(opts)
if ('-h' in dopt) or ('--help' in dopt):
	print usage()
	sys.exit(0)
	
nflnfingtchain = dopt['-G']
dirout = dopt['-o']
chain1ext = dopt.get('--chain.ext', 'mb.run1.t')
nbchains = int(dopt.get('--nb.chains', 2))
nbthreads = int(dopt.get('--threads', 1))
verbose = (('--verbose' in dopt) or ('-v' in dopt))

with open(nflnfingtchain, 'r') as flnfingtchain:
	lnfingtchain = [line.strip('\n') for line in flnfingtchain]

def parseChainSetArgs(nfingtchain1):
	bngt, extgt = os.path.basename(nfingtchain1).split('.', 1)
	fam = bngt.rsplit('-', 1)[0]
	nfoutcolGtrees = "%s/%s-Gtrees.nwk"%(dirout, bngt)
	nfingtchains = [nfingtchain1.replace(chain1ext, chain1ext.replace('1', str(k))) for k in range(1, nbchains+1)]
	parseChain(nfingtchains, nfchainout=nfoutcolGtrees, verbose=verbose)


# run in sequential
if nbthreads==1:
	for nfingtchain in lnfingtchain:
		if verbose: print nfingtchain
		parseChainSetArgs(nfingtchain)
# or run in parallel
else:
	pool = mp.Pool(processes=nbthreads)
	res = pool.map(parseChainSetArgs, lnfingtchain)
