#!/usr/bin/env python2.7

import xml.etree.ElementTree as ET

import sys, os, getopt

sys.setrecursionlimit(20000)

deventshort = {'speciation':'S', 'leaf':'S', 'duplication':'D',  'loss':'L',  'transfer':'T'}

def parse_element(elt, fatelt=None):
	levt = []
	# record who was the father
	elt.father = fatelt
	# use pre-order traversal
	if 'eventsRec' in elt.tag:
		evt = parse_eventsrec_element(elt)
		if evt[0] == 'branchingOut':
			assert ('clade' in fatelt.tag)
			fatelt.bo = evt[1]
		elif evt[0] == 'transferBack':
			# look for a branchingOut attribute on the grand father element
			ori = elt.father.father.bo
			assert (ori is not None)
			levt.append(('T', evt[1], ori))
		else:
			levt.append((deventshort[evt[0]], evt[1]))
	for child in elt:
		levt += parse_element(child, fatelt=elt)
	return levt
	
def parse_eventsrec_element(elt):
	evtelt = elt[0]
#	print evtelt.__dict__
	try:
		eloc = evtelt.attrib['speciesLocation']
	except KeyError:
		eloc = evtelt.attrib['destinationSpecies']
	etype = evtelt.tag.rsplit('}', 1)[1]
	return (etype, eloc)

nfinxml = sys.argv[1]
nfout = sys.argv[2]

parser = ET.XMLParser(encoding="utf-8")
elttree = ET.parse(nfinxml, parser=parser)
elttreeroot = elttree.getroot()

recgt = None
for child in elttreeroot:
	if 'recGeneTree' in child.tag:
		recgt = child
		break

levt = parse_element(recgt)
with open(nfout, 'w') as fout:
	for evt in levt:
		fout.write('\t'.join(evt)+'\n')
