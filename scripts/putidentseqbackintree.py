#!/usr/bin/python
#~ import os
import sys, tree2, getopt

sys.setrecursionlimit(20000)

def checkBS(tree, maxNone=2):
	"""count how many node ave no branch support documented; raise an error when above threshold
	
	by default allow a maximum of 2 without support to account for branches created when rooting
	"""
	nnodenobs = 0
	for node in tree:
		if node.is_leaf(): continue
		if node.bs() is None: nnodenobs += 1
	
	if nnodenobs > maxNone:
		raise ValueError, "too many (%d) nodes without branch support documented"%nnodenobs


opts, args = getopt.getopt(sys.argv[1:], 'h', ['input.nr.tree=', 'output.tree=', 'list.identical.seq=', 'ref.rooted.nr.tree=', 'help'])

dopt = dict(opts)

if ('-h' in dopt) or ('--help' in dopt):
	#~ print usage()
	exit()

nfintree = dopt['--input.nr.tree']
nfidentseqs = dopt.get('--list.identical.seq')
nfrefroottree = dopt.get('--ref.rooted.nr.tree')
nfouttree = dopt.get('--output.tree', nfintree+'.full')


intree = tree2.Node(file=nfintree)
try:
	checkBS(intree)
except ValueError:
	print "could not find branch supports, look in comments field"
	intree = tree2.Node(file=nfintree, keep_comments=True)
	for n in intree:
		if str(n.comment()).isdigit():
			n.set_bs(float(n.comment()))
	checkBS(intree)

dseq = {}
if nfidentseqs:
	with open(nfidentseqs, 'r') as fidentseqs:
		for line in fidentseqs:
			lsp = line.rstrip('\n').split('\t')
			dseq.setdefault(lsp[0], []).append(lsp[1])

if nfrefroottree:
	refroottree = tree2.Node(file=nfrefroottree)
	subrootchildren = refroottree.get_children()
	# map sub-root clades in reference tree to the input tree
	mappedclades = [intree.map_to_node(src.get_leaf_labels()) for src in subrootchildren]
	# find outgroup index as the one not being root
	for i, mc in enumerate(mappedclades):
		if not mc.is_root(): break
	outgroup = mappedclades[i]
	# branch length under the root such as: [bl(outgroup), bl(other clade)]
	subrootbrlens = [subrootchildren[t].lg() for  t in (i, 1-i)]
	# re-root the input tree as in the reference
	intree.newOutgroup(outgroup, branch_lengths=subrootbrlens)


for nrseq in dseq:
	leaf = intree[nrseq]
	if not leaf: raise IndexError, 'the nr sequence named %s is missing from the tree!'%nrseq
	newfat = leaf.create_node_above(newboot=100.0)
	for identseq in dseq[nrseq]:
		newfat.create_leafnode_below(label=identseq)
	
	
intree.write_newick(nfouttree)
