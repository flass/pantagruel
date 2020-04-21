#!/usr/bin/env python2.7

import tree2

import sys, os, getopt

#sys.setrecursionlimit(20000)

orinhxcommsep=':'
cleannhxcommsep='/'

def parse_branch_annot(annot):
	levt = []
	evloc = None
	for elt in annot.split(cleannhxcommsep):
		if elt=='&&NHX': continue
		e1, e2 = elt.split('=')
		if e1=='S':
			evloc = e2
		else:
			if e2!='N':
				if e1=='D':
					levt.append(('D', evloc, ''))
				elif e1=='H':
					don, rec = e2.split('@')[1:]
					levt.append(('T', rec, don))
	if not levt:
		# no D or H/T event recorded: just a speciation S
		levt.append(('S', evloc, ''))
	return levt
	


def parseMultiNHXRec(nfinnhx):
	lrecgt = []
	devtcount = {}
	with open(nfinnhx, 'r') as finnhx:
		# first clean this dirty NHX format with semilocons within comment that messes the parsing
		# replace in-comment ':' with '/'
		for line in finnhx:
			scleannhx = ''
			incomment = 0
			for char in line:
				if char=='[': incomment += 1
				elif char==']': incomment -= 1
				if incomment:
					if char==orinhxcommsep: char=cleannhxcommsep
				scleannhx += char

	#		print scleannhx
			recgt = tree2.AnnotatedNode(newick=scleannhx, keep_comments=True)
			recgt.complete_node_ids()
			lrecgt.append(recgt)
			for node in recgt:
				for evt in parse_branch_annot(node.comment()):
					devtcount[evt] = devtcount.setdefault(evt, 0) + 1
	return (devtcount, lrecgt)

if __name__ == "__main__":

	nfinnhx = sys.argv[1]
	nfout = sys.argv[2]

	devtcount, lgt = parseMultiNHXRec(nfinnhx)

	with open(nfout, 'w') as fout:
		for evt, n in devtcount.iteritems():
			fout.write('\t'.join(evt+(str(n),))+'\n')
