#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""makes a matrix of collapsed clade (CC) occurence in genome populations ; only to use if not running Pantagruel task 06.5"""

import os, sys, getopt
import multiprocessing as mp

def parse_pops(nfpopdef):
	dcode2pop = {}
	spop = set([])
	with open(nfpopdef, 'r') as fpopdef:
		for line in fpopdef:
			if line.startswith('#'): continue
			lsp = line.rstrip('\n').split('\t')
			p = lsp[0]
			spop.add(p)
			for code in lsp[1:]:
				dcode2pop[code] = p
	lpop = sorted(spop)
	return (lpop, dcode2pop)

def parse_fam_constraint_folder(dfamcons):
	fam = os.path.basename(dfamcons)
	lnfcons = os.listdir(dfamcons)
	dmatCCcount = {}
	for nfcons in lnfcons:
		with open(os.path.join(dfamcons, nfcons)) as fcons:
			for line in fcons:
				if not line.startswith('constraint'): continue
				cladef, cdscodes = line.rstrip('\n').split(' = ')
				famcla = fam+'-'+cladef.split()[1]
				lcdscodes = cdscodes.split()
				dmatCCrow = {}
				for cdscode in lcdscodes:
					p = cdscode.split('_',1)[0]
					dmatCCrow[p] = dmatCCrow.get(p, 0) + 1
				dmatCCcount[famcla] = dmatCCrow
	return dmatCCcount

def main():
	opts, args = getopt.getopt(sys.argv[1:], 'c:p:o:hvT:', ['CC_defs=', 'populations=', 'out=', 'threads=', 'help', 'verbose']) #, 'reuse=', 'max.recursion.limit=', 'logfile='
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	ccdefdir = dopt.get('-c', dopt.get('--CC_defs'))
	nfpopdef = dopt.get('-p', dopt.get('--populations'))
	nfmatout = dopt.get('-o', dopt.get('--out'))
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	nbthreads = int(dopt.get('--threads', dopt.get('-T', -1)))
	if nbthreads < 1:
		nbthreads = mp.cpu_count()
	
	if verbose:
		print "# call: %s"%(' '.join(sys.argv))
		print "# dopt:", dopt
	
	# get population definitions
	lpop, dcode2pop = parse_pops(nfpopdef)
	
	iterdfamcons = (os.path.join(ccdefdir, dfamcons) for dfamcons in os.listdir(ccdefdir))
	
	# parse the collapsed clade definitions in parallel
	if nbthreads == 1:
		ildevents = (parse_fam_constraint_folder(dfamcons) for dfamcons in iterdfamcons)
	else:
		pool = mp.Pool(processes=nbthreads)
		# an iterator is returned by imap(); one needs to actually iterate over it to have the pool of parrallel workers to compute
		idmatCCcount = pool.imap_unordered(parse_fam_constraint_folder, iterdfamcons)
	
	fmatout = open(nfmatout, 'w')
	for dmatCCcount in idmatCCcount:
		for famcla, dmatCCrow in dmatCCcount.iteritems():
			fmatout.write('\t'.join([famcla]+[str(dmatCCrow.get(p, 0)) for p in lpop]))
	
	fmatout.close()

def usage():
	s = "Usage:\n"
	s += "python %s --CC_defs _mb_constraint_folder_ --out _output_matrix_file_\n"%sys.argv[0]
	s += "Options:\n"
	s += "  -c|--CC_defs      path_to_folder   path to input folder containing a folder per gene family,\n"
	s += "                                       each containing a MrBayes constraint file per collapsed clade\n"
	s += "  -p|--populations  path_to_file     path to input file defining the populations\n" 
	s += "                                       (a tab-separated table, with left field the population name\n"
	s += "                                       and right field the space-separated list of genome codes)\n"
	s += "  -o|--CC_defs      path_to_file     path to output matrix file\n"
	return s

################## Main execution

if __name__=='__main__':
	
	main()