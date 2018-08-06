#!/usr/bin/python
import sys, os
import tree2
from parseALErec import parseALERecFile
import copy

scaleFreqToWidth = 500

presentColor = (200, 200, 255)
transferColor = (200, 200, 80)#, 127)

recfilesuffix = "-collapsed.ale.uml_rec"

def linesplit(line):
	return [field.strip('" ') for field in line.rstrip('\n').split('\t')]

if len(sys.argv)<3:
	print "Usage: %s /path/to/lineage_module_event_table /path/to/reference_tree /path/to/output_folder [/path/to/reconciliation_folder]"
	sys.exit(2)

nflnflineagecommevents = sys.argv[1]
nfrefspetree = sys.argv[2]
dirout = sys.argv[3]
# optional argument to give context of gene family occurrence to gene lineage
if len(sys.argv)>4:
	dirrec = sys.argv[4]
else:
	dirrec = None

refspetree = tree2.AnnotatedNode(file=nfrefspetree)
			
with open(nflnflineagecommevents, 'r') as flnflineagecommevents:
	lnflineagecommevents = [line.rstrip('\n') for line in flnflineagecommevents]

dfamspetree = {}
for nflineagecommevents in lnflineagecommevents:
	flineagecommevents = open(nflineagecommevents, 'r')
	lineagecomm = os.path.basename(nflineagecommevents).rsplit('.', 1)[0]
	dirlineageout = os.path.join(dirout, lineagecomm)
	if not os.path.isdir(dirlineageout): os.mkdir(dirlineageout)
	curfamily = None
	curlineage = None
	curspetree = None
	dnodefreq = {}
	ltrans = []
	header = linesplit(flineagecommevents.readline())
	for line in flineagecommevents:
		family, lineage, event, freq, evtype, reclabel, donlabel = linesplit(line)
		if family!=curfamily:
			if dirrec:
				# load pobability density of gene presence of the whole gene family over the species tree
				# !!! when reconciliation used partially collapsed species tree, requires a matching of uncollapsed to collapsed nodes (NOT IMPLEMENTED)
				if family not in dfamspetree:
					nfrec = os.path.join(dirrec, "%s%s"%(family, recfilesuffix))
					recspetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt = parseALERecFile(nfrec, reftreelen=refspetree)
					for node in recspetree:
						node.branchwidth = float(dnodeevt[node.label()][-1])/scaleFreqToWidth
					dfamspetree[family] = recspetree
				else:
					recspetree = dfamspetree[family]
			curfamily = family
		if lineage!=curlineage:
			if not (curlineage is None):
				# write out previous lineage projection
				nfoutrad = os.path.join(dirlineageout, "lineage_%s_projection"%curlineage)
				curspetree.write_newick(nfoutrad+".nwk", ignoreBS=False)
				curspetree.writeSvgTree(nfoutrad+".svg", padleaves=True, supports=False, phylofact=10000, branchwidths='branchwidth', \
		                                  treetype='species', transfers=ltrans, textorbit=5, modstyle="stroke-width:1; ", \
		                                  transfercolor=transferColor, transferpathtype='arc', transferwidth='freq')
			if dirrec: curspetree = copy.deepcopy(recspetree)
			else: curspetree = copy.deepcopy(refspetree)
			for node in curspetree:
				node.set_bs(0.0)
			dnodefreq = {}
			ltrans = []
			curlineage = lineage
		
		scaledfreq = float(freq)/scaleFreqToWidth
		if evtype=='T':
			# transfer events are incorectly soted in PanteroDB v0.3: donor and recipient labels are swapped;
			# here is a quick fix for graphic representation, not intended to last!!
			refrec = curspetree[donlabel]
			refdon = curspetree[reclabel]
			ltrans.append((reclabel, donlabel, scaledfreq))
			#~ ltrans.append((donlabel, reclabel, scaledfreq))
		else:
			refrec = curspetree[reclabel]
			refdon = None
		refrec.set_bs(refrec.bs()+scaledfreq)
		refrec.branchwidth = scaledfreq
		refrec.edit_color(presentColor)
	
	flineagecommevents.close()
