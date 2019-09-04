#!/usr/bin/python
#~ import os
import sys, tree2, getopt

sys.setrecursionlimit(20000)

opts, args = getopt.getopt(sys.argv[1:], 'h', ['input.nr.tree=', 'output.tree=', 'list.identical.seq=', 'ref.rooted.nr.tree=', 'help'])

dopt = dict(opts)

if ('-h' in dopt) or ('--help' in dopt):
	#~ print usage()
	exit()

nfintree = dopt['--input.nr.tree']
nfidentseqs = dopt.get('--list.identical.seq')
nfrefroottree = dopt.get('--ref.rooted.nr.tree')
nfouttree = dopt.get('--output.tree', nfintree+'.full')

intree = tree2.read_check_newick(nfintree)

dseq = {}
if nfidentseqs:
	with open(nfidentseqs, 'r') as fidentseqs:
		for line in fidentseqs:
			lsp = line.rstrip('\n').split('\t')
			dseq.setdefault(lsp[0], []).append(lsp[1])

if nfrefroottree:
	refroottree = tree2.read_check_newick(nfrefroottree)
	
	if not refroottree.is_bifurcated():
		raise ValueError, "the provided 'rooted' tree is not bifurcated"
	nnodenobs = tree2.checkBS(refroottree, maxNoBS=-1)
	if nnodenobs <4:
		# the rooted tree can be used as is
		intree = refroottree
	else:
		print "must transfer supports from the tree with bipartitions"
		subrootchildren = refroottree.get_children()
		# map sub-root clades in reference tree to the input tree
		mappedclades = [intree.map_to_node(src.get_leaf_labels()) for src in subrootchildren]
		# find outgroup index as the one not being root
		for i, mc in enumerate(mappedclades):
			if not mc.is_root():
				break
		else:
			raise IndexError, "Could not find in '%s' the basal clade forming an outgroup in '%s'"%(nfintree, nfrefroottree)
		outgroup = mappedclades[i]
		outfat = outgroup.go_father()

		# branch length under the root such as: [bl(outgroup), bl(other clade)]
		subrootbrlens = [subrootchildren[t].lg() for  t in (i, 1-i)]
		# re-root the input tree as in the reference
		intree.newOutgroup(outgroup, branch_lengths=subrootbrlens)
		intree.resolveNode([outgroup])
		intree = intree.go_root()
		assert intree.is_bifurcated()
		if subrootchildren[0].lg() > 0:
			# correct the branch lengths
			for srnode in intree.get_children():
				for src in subrootchildren:
					# find the matching reference node
					if set(src.get_leaf_labels()) ==  set(srnode.get_leaf_labels()):
						srnode.set_lg(src.lg())
						break
				else:
					raise IndexError, "could not mnatch the reference in input tree nodes"

for nrseq in dseq:
	leaf = intree[nrseq]
	if not leaf: raise IndexError, 'the nr sequence named %s is missing from the tree!'%nrseq
	newfat = leaf.create_node_above(newboot=100.0)
	for identseq in dseq[nrseq]:
		newfat.create_leafnode_below(label=identseq)
	
	
intree.write_newick(nfouttree)
