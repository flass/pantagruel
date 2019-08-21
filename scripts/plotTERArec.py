#!/usr/bin/python
import tree2
import sys, os, getopt
from parseTERArec import *

def main(nfrec, nfreftree, maxrecgt=1, recformat='tera', sgsep='_', phylofact=1000.0, restrictclade=None, verbose=False):
	reftree = tree2.AnnotatedNode(file=nfreftree, namesAsNum=True)
	i = -1	
	for dnodefreq, dlevt in parseTERARecFile(nfrec, recformat=recformat, sgsep=sgsep, verbose=verbose):
		i += 1
		# write SVG species tree
		nfoutspe = '%s_%d_maprec2spetree.svg'%(nfrec, i)
		if restrictclade: st = reftree.restrictToLeaves(restrictclade)
		else: st = reftree
		print dnodefreq
		lleaffreq = [(lab, f) for lab, f in dnodefreq.items() if st[lab].is_leaf()]
		st.writeSvgTree(nfoutspe, padleaves=True, supports=False, phylofact=phylofact, branchwidths=dnodefreq, textorbit=5, \
		 treetype='species', transfers=dlevt['T'], duplications=dlevt['D'], losses=dlevt['L'], counts=lleaffreq, \
		 transferwidth='freq', modstyle="stroke-width:1; ", padstyle="stroke:red; stroke-width:0.5; stroke-dasharray:1,1; ")
		# transfercolor='green', 
		print os.path.basename(nfoutspe)

def usage():
	s  = 'python plotTERArec.py [options] reconcilationFile{.mr|.txt} [reconcilation2.ale.uml_rec [... reconcilationN.ale.uml_rec]]\n'
	s += 'Mandatory argument:\n'
	s += '\t--reftree path\tpath to tree file in newick format; tree must be the same than that was used in the TERA reconciliation\n'
	s += 'Options:\n'
	s += '\t\t allows to have branch lengths on species tree plot.\n'
	s += '\t--maxrecgt int\tmaximum number of alternative reconcilaition scenarios (among the top ML ones) to generate a plot for\n'
	s += '\t\t defaults to 1.\n'
	s += '\t--species-gene-separator str\tstring separator used on gene tree leaves to separate species tag from gene id; default: \'_\'.\n'
	#~ s += '\t--restrict-to-clade str\tstring[,string[,...]] label(s) of node(s) of the species tree to which the plots will be restricted; default: none.\n'
	s += '\t---magnify scaling factor for the tree graphics\n'
	s += '\t---format the reconcilaition file format: either \'mowgli\' or \'tera\'; files are normally coming out of ecceTERA\n'
	s += '\t          with .mr and .txt extensions, respectively; default: guess form file extension, or raise an error if not possible.\n'
	s += '\t-h | --help print this help message and exit.\n'
	return s

## execute script
if __name__ == "__main__":
	## arg parsing
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hv", ['reftree=', 'maxrecgt=', 'species-gene-separator=', 'magnify=', 'format=', 'help', 'verbose']) #, 'restrict-to-clade='
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err)  # will print something like "option -a not recognized"
		print usage()
		sys.exit(2)
		
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	nfreftree = dopt['--reftree']
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	
	maxrecgt = int(dopt.get('--maxrecgt', 1))
	recformat = dopt.get('--format')
		
	sgsep = dopt.get('--species-gene-separator', '_')
	phylofact = float(dopt.get('--magnify', 1000))
	restrictclade = dopt.get('--restrict-to-clade')
		
	lnfrec = args
	if len(lnfrec)<1: raise ValueError, "need at least one argument (file path[s])"
	
	for nfrec in lnfrec:
		if not recformat:
			if nfrec.endswith('.mr'): recformat = 'mowgli'
			elif nfrec.endswith('.txt'): recformat = 'tera'
			else: raise ValueError, "could not guess the format of the file from its extension; please specify either 'mowgli' or 'tera'."
		main(nfrec, nfreftree, maxrecgt, recformat, sgsep, phylofact, restrictclade, verbose=verbose)
