#!/usr/bin/python
from parseALErec import getOrthologues
import ptg_utils as ptg
import tree2
import os, glob, sys, getopt

methods = ['strict', 'unicopy', 'mixed']

def summaryOGs(ogs, dlabs, N):
	n = len(ogs)
	print "number of OGs:", n
	cov = sum([len(x) for x in ogs])
	print "coverage of leaves:", cov, '/', N
	if cov != N:
		raise ValueError, "unclassified leaves:%s"%(repr(set(dlabs.values()) - set(reduce(lambda x,y: list(x)+list(y), ogs, []))))
	return str(n)

def orthoFromSampleRecs(nfrec, foutdiffog=None, nsample=[], methods=['mixed'], colourTree=False, **kw):
	""""""
	verbose = kw.get('verbose')
	fam = os.path.basename(nfrec).split('-', 1)[0]
	if verbose: print "\n# # # %s"%fam
	# collect the desired sample from the reconciliation file
	dparserec = parseALERecFile(nfrec, skipLines=True, skipEventFreq=True, nsample=nsample, returnDict=True)
	lrecgt = dparserec['lrecgt']
	if userefspetree:
		refspetree = dparserec['spetree']
	else:
		refspetree = None
	
	ddogs = {}
	dlabs = {}
	dnexuslabev2num = {}
	dnexustrans = {}
	for i, recgenetree in enumerate(lrecgt):
		if nsample: g = nsample[i]
		else: g = i
		if verbose: print recgenetree
		if verbose: print "\n# # reconciliation sample %d"%g
		N = recgenetree.nb_leaves()
		if set(['strict', 'mixed']) & set(methods):
			if verbose: print "\n# strict_ogs:\n"
			strict_ogs, dlabs = getOrthologues(recgenetree, method='strict', refspetree=refspetree, dlabs=dlabs, **kw)
			n1 = str(summaryOGs(strict_ogs, dlabs1, N))
		else:
			strict_ogs = None, n1 = 'NA'
		if 'unicopy' in methods:
			if verbose: print "\n# unicopy_ogs:\n"
			unicopy_ogs, dlabs = getOrthologues(recgenetree, method='unicopy', refspetree=refspetree, dlabs=dlabs, **kw)
			n2 = summaryOGs(unicopy_ogs, dlabs2, N)
		else:
			unicopy_ogs = None, n2 = 'NA'
		if 'mixed' in methods:
			if verbose: print "\n# mixed_ogs:\n"
			mixed_ogs, dlabs = getOrthologues(recgenetree, method='mixed', gain_ogs=strict_ogs, refspetree=refspetree, dlabs=dlabs, **kw) #
			n3 = summaryOGs(mixed_ogs, dlabs3, N)
		else:
			mixed_ogs = None, n3 = 'NA'
		
		if foutdiffog or verbose: 
			o12 = str(sum([int(o in strict_ogs) for o in unicopy_ogs])) if (strict_ogs and unicopy_ogs) else 'NA'
			o13 = str(sum([int(o in strict_ogs) for o in mixed_ogs])) if (strict_ogs and mixed_ogs) else 'NA'
			o23 = str(sum([int(o in unicopy_ogs) for o in mixed_ogs])) if (mixed_ogs and unicopy_ogs) else 'NA'
		if verbose:
			print "\n# summary:\n"
			print "overlap strict_ogs with unicopy_ogs:", o12
			print "overlap strict_ogs with mixed_ogs:", o13
			print "overlap unicopy_ogs with mixed_ogs:", o23
		if foutdiffog:
			foutdiffog.write('\t'.join([fam, g, n1, n2, n3, o12, o13, o23])+'\n')
		
		if colourTree:
			if not dnexuslabev2num:
				for l, leaf in enumerate(recgenetree.get_leaves()):
					labev = leaf.label()
					lab = dlabs[labev]
					dnexuslabev2num[labev] = l
					dnexustrans[l] = lab
				ltaxnexus = [dnexustrans[l] for l in range(len(dnexustrans))]
			for node in recgenetree:
				if node.is_leaf():
					node.edit_label(str(dnexuslabev2num[node.label()]))
				else:
					node.edit_label('')
		ddogs[g] = {'strict':strict_ogs, 'unicopy':unicopy_ogs, 'mixed':mixed_ogs, 'recgenetree':recgenetree}
		if verbose: print "\n# # # # # # # #"
	
	gs = namsample if namsample else range(len(lrecgt))
	for method in methods:
		ltrees = []
		nfoutrad = os.path.join(outortdir, method, "%s_%s"%(fam, method))
		if colourTree:
			for g in gs:
				recgenetree = ddogs[g]['recgenetree']
				ogs = ddogs[g][method]
				ptg.colour_tree_with_leaf_groups(recgenetree, ogs)
				#~ ptg.colour_tree_with_constrained_clades(recgenetree, ogs)
				ltrees.append(recgenetree)
			tree2.write_nexus(ltrees, nfout=nfoutrad+"_orthologous_groups.nex", ltax=ltaxnexus, dtranslate=dnexustrans, ignoreBS=True)
		with open(nfoutrad+".orthologs", 'w') as foutort:
			for g in gs:
				ogs = ddogs[g][method]
				foutort.write(';'.join([','.join(x) for x in ogs])+'\n')

def main(alerecdir, outortdir, summaryOverlap=False, **kw):
	
	verbose = kw.get('verbose')
	methods = kw.get('methods', ['mixed'])
	lnfrec = glob.glob('%s/*ale.ml_rec'%(alerecdir))

	for d in methods:
		pd = os.path.join(outortdir, d)
		if not os.path.isdir(pd):
			os.mkdir(pd)
	
	if summaryOverlap:
		foutdiffog = open(os.path.join(outortdir, 'diff_ortho_methods'), 'w')
		foutdiffog.write('\t'.join(['family', 'sampled_rec', 'nOG_strict', 'nOG_unicopy', 'nOG_mixed', 'overlap_strict_unicopy', 'overlap_strict_mixed', 'overlap_unicopy_mixed'])+'\n')
	else:
		foutdiffog = None
	
	for nfrec in lnfrec:
		orthoFromSampleRecs(nfrec, foutdiffog=foutdiffog, **kw)
		if verbose: print " # # # # # # # \n"

	if summaryOverlap: foutdiffog.close()

def usage():
	s =  'python %s -i /path/to/input.reconciliation.dir -o /path/to/output.dir [OPTIONS]\n'%sys.argv[0]
	s += ' Options:\n'
	s += '  --verbose {0,1,2}\tverbose mode, from none to plenty.\n'
	s += '  -v\t\t\tequivalent to --verbose=1.\n'
	return s

if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:hv', ['input.reconciliation.dir=', 'output.dir=', 'methods=', \
														'summary', 'max.frac.extra.spe=', \
														'threads=', 'verbose=', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	alerecdir = dopt.get('--input.reconciliation.dir', dopt.get('-i'))
	outortdir = dopt.get('--output.dir', dopt.get('-o'))
	if not (alerecdir and outortdir):
		print "missing mandatory arguments"
		print usage()
		sys.exit(2)
	
	userefspetree = dopt.get('--use.species.tree')
	verbose = int(dopt.get('--verbose', ('-v' in dopt)))
	ALEmodel = dopt.get('--ALE.model', 'dated')
	summaryOverlap = ('--summary' in dopt)
	trheshExtraSpe = float(dopt.get('--max.frac.extra.spe', 0.1))
	
	main(alerecdir, outortdir, summaryOverlap=summaryOverlap, ALEmodel=ALEmodel, userefspetree=userefspetree, \
	     trheshExtraSpe=trheshExtraSpe, verbose=verbose)
