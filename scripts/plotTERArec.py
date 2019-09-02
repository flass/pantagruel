#!/usr/bin/python
import tree2
import sys, os, getopt
from parseTERArec import *

def main(nfrec, nfreftree, nfgenetree, maxrecgt=1, recformat='tera', sgsep='_', phylofact=1000.0, restrictclade=None, verbose=False, **kw):
	try:
		genetree = tree2.Node(file=nfgenetree, namesAsNum=True)
	except ValueError:
		genetree = tree2.Node(file=nfgenetree, namesAsNum=True, branch_lengths=False)
	reftree = tree2.AnnotatedNode(file=nfreftree, namesAsNum=True)
	if restrictclade: st = reftree.restrictToLeaves(restrictclade)
	else: st = reftree
	# check presence of outgroup/dead lineage branch if necessary
	if recformat=='tera':
		if not (kw.get('noDeadStories') or (deadlabnum in st.get_leaf_labels())):
			if (outtaxlab in st.get_leaf_labels()):
				# must adapt mowgli-compliant species tree
				st[outtaxlab].edit_label(deadlabnum)
			else:
				maxd = reftree.max_leaf_distance()
				outgroup = tree2.AnnotatedNode(lleaves=[deadlabnum])
				outgroup.get_children()[0].set_lg(maxd*3)
				outgroup.link_child(reftree, newlen=maxd*2)
				reftree = outgroup
				reftree.complete_internal_labels(prefix='')
#			else:
#				raise ValueError, "the provided species tree should feature a branch labaelled 'OUTGROUP' or '-1' to represent the dead/unsampled lineages"
	elif recformat=='mowgli':
		if not (outtaxlab in st.get_leaf_labels()):
			if (deadlabnum in st.get_leaf_labels()):
				# must adapt mowgli-compliant species tree
				st[deadlabnum].edit_label(outtaxlab)
			else:
				outgroup = tree2.AnnotatedNode(lleaves=[outtaxlab])
				outgroup.get_children()[0].set_lg(maxd*3)
				outgroup.link_child(reftree, newlen=maxd*2)
				reftree = outgroup
				reftree.complete_internal_labels(prefix='')
#			else:
#				raise ValueError, "the provided species tree should feature a branch labaelled 'OUTGROUP' or '-1' to represent the dead/unsampled lineages"
	for i, rec in enumerate(parseTERARecFile(nfrec, genetree=genetree, recformat=recformat, sgsep=sgsep, verbose=verbose, **kw)):
		dnodefreq, dlevt = rec
		# write SVG species tree
		tag = '_no_dead' if kw.get('noDeadStories') else ''
		nfoutspe = '%s_%d_maprec2spetree%s.svg'%(nfrec, i, tag)
		lleaffreq = [(lab, f) for lab, f in dnodefreq.items() if st[lab].is_leaf()]
		st.writeSvgTree(nfoutspe, padleaves=True, supports=False, phylofact=phylofact, branchwidths=dnodefreq, textorbit=5, \
		 treetype='species', transfers=dlevt['T'], duplications=dlevt['D'], losses=dlevt['L'], counts=lleaffreq, \
		 transferwidth='freq', modstyle="stroke-width:1; ", padstyle="stroke:red; stroke-width:0.5; stroke-dasharray:1,1; ")
		# transfercolor='green', 
		print os.path.basename(nfoutspe)

def usage():
	s  = 'python plotTERArec.py [options] reconcilationFile{.mr|.txt} [reconcilation2.ale.uml_rec [... reconcilationN.ale.uml_rec]]\n'
	s += 'Mandatory argument:\n'
	s += '\t--reftree or --speciestree path\tpath to S tree file in newick format\n'
	s += '\t--genetree path\tpath to G tree file in newick format\n'
	s += '\t\tS and G trees must be the same than those that were used in the TERA reconciliation, or equivalent to them (e.g. as provided in ecceTERA output)\n'
	s += 'Options:\n'
	s += '\t--no-dead-stories\tsimplify the sceanrio by directly connecting \'live\' lineages of the species tre, ignoring events that took place among the \'dead\'.\n'
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
		opts, args = getopt.getopt(sys.argv[1:], "hv", ['reftree=', 'speciestree=', 'genetree=', 'maxrecgt=', 'species-gene-separator=', 'magnify=', 'format=', 'no-dead-stories', 'help', 'verbose']) #, 'restrict-to-clade='
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err)  # will print something like "option -a not recognized"
		print usage()
		sys.exit(2)
		
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	nfreftree = dopt.get('--reftree', dopt.get('--speciestree'))
	if not nfreftree: raise ValueError, "reference tree file argument is missing; please specify it with either --reftree or --speciestree"
	nfgenetree = dopt['--genetree']
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	
	maxrecgt = int(dopt.get('--maxrecgt', 1))
	recformat = dopt.get('--format')
		
	sgsep = dopt.get('--species-gene-separator', '_')
	phylofact = float(dopt.get('--magnify', 1000))
	restrictclade = dopt.get('--restrict-to-clade')
	noDeadStories = ('--no-dead-stories' in dopt)
		
	lnfrec = args
	if len(lnfrec)<1: raise ValueError, "need at least one argument (file path[s])"
	
	for nfrec in lnfrec:
		if not recformat:
			if nfrec.endswith('.mr'): recformat = 'mowgli'
			elif nfrec.endswith('.txt'): recformat = 'tera'
			else: raise ValueError, "could not guess the format of the file from its extension; please specify either 'mowgli' or 'tera'."
		main(nfrec, nfreftree, nfgenetree, maxrecgt=maxrecgt, recformat=recformat, sgsep=sgsep, phylofact=phylofact, restrictclade=restrictclade, noDeadStories=noDeadStories, verbose=verbose)
