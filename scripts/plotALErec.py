#!/usr/bin/python
import tree2
import sys, os, getopt
import re
from parseALErec import *
#~ import cPickle

def main(nfrec, reftreelen=None, maxrecgt=1, sgsep='_', restrictclade=None):
	spetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt = parseALERecFile(nfrec, reftreelen=reftreelen, restrictclade=restrictclade)
	dnodefreq = dict([(node.label(), float(dnodeevt[node.label()][-1])) for node in spetree])
	nsample = len(lrecgt)
	recgtsample = ''.join(recgtlines)
	# parse reconciled gene trees
	# and extract (exact) event-wise event frequency
	dexactevt = {}
	for i in range(min(len(lrecgt), maxrecgt)):
		recgt = lrecgt[i]
		# gather scenario-scpecific events
		dlevt, dnodeallevt = parseRecGeneTree(recgt, spetree, dexactevt, recgtsample, nsample, sgsep=sgsep, restrictlabs=restrictlabs)
		# gather copy number
		dcopynum = {}
		for leaflab in recgt.get_leaf_labels():
			spe = leaflab.split(sgsep)[0]
			if restrictlabs and not (spe in restrictlabs): continue
			dcopynum[spe] = dcopynum.setdefault(spe, 0) + 1
		lcopynum = dcopynum.items()
		# write SVG species tree
		nfoutspe = '%s_%d_maprec2spetree.svg'%(nfrec, i)
		if restrictclade: st = subspetree
		else: st = spetree
		st.writeSvgTree(nfoutspe, padleaves=True, supports=False, phylofact=(10000 if reftreelen else 10), branchwidths=dnodefreq, textorbit=5, \
		 treetype='species', transfers=dlevt['T'], duplications=dlevt['D'], losses=dlevt['L'], counts=lcopynum, \
		 modstyle="stroke-width:1; ", padstyle="stroke:red; stroke-width:0.5; stroke-dasharray:1,1; ")
		print os.path.basename(nfoutspe)
		# write SVG species tree
		nfoutrec = '%s_%d_recgenetree.svg'%(nfrec, i)
		#~ if restrictlabs: gt = recgt.map_to_node(restrictlabs, force=True, useSpeDict=True)
		#~ else: gt = recgt
		recgt.writeSvgTree(nfoutrec, padleaves=False, phylofact=1000, textorbit=10, \
		 treetype='reconciliation', nodeevents=dnodeallevt, \
		 modstyle="stroke-width:1; ", dxindex='nodeid')
		 # losses=dnodeallevt['L'] \
		print os.path.basename(nfoutrec)

def usage():
	s  = 'python plotALErec.py [options] reconcilation1.ale.uml_rec [reconcilation2.ale.uml_rec [... reconcilationN.ale.uml_rec]]\n'
	s += 'Options:\n'
	s += '\t--reftree path\tpath to tree file in newick format; tree must be the same than that used in the ALE reconciliation\n'
	s += '\t\t allows to have branch lengths on species tree plot.\n'
	s += '\t--maxrecgt int\tmaximum number of alternative reconcilaition scenarios (among the top ML ones) to generate a plot for\n'
	s += '\t\t defaults to 1.\n'
	s += '\t--species-gene-separator str\tstring separator used on gene tree leaves to separate species tag from gene id; default: \'_\'.\n'
	s += '\t--restrict-to-clade str\tstring[,string[,...]] label(s) of node(s) of the species tree to which the plots will be restricted; default: none.\n'
	return s

## execute script
if __name__ == "__main__":
	## arg parsing
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h", ['reftree=', 'maxrecgt=', 'species-gene-separator=', 'restrict-to-clade=', 'help'])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err)  # will print something like "option -a not recognized"
		print usage()
		sys.exit(2)
		
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	if '--reftree' in dopt:
		nfreftreelen = dopt['--reftree']
		reftreelen = tree2.Node(fic=nfreftreelen)
		print "Using branch lengths from reference tree '%s'"%(nfreftreelen)
	else:
		reftreelen = None
	if '--maxrecgt' in dopt: maxrecgt = int(dopt['--maxrecgt'])
	else: maxrecgt = 1
	sgsep = dopt.get('--species-gene-separator', '_')

	restrictclade = dopt.get('--restrict-to-clade')
		
	lnfrec = args
	if len(lnfrec)<1: raise ValueError, "need at least one argument (file path[s])"
	
	for nfrec in lnfrec:
		main(nfrec, reftreelen, maxrecgt, sgsep, restrictclade)
