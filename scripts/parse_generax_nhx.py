#!/usr/bin/env python2.7

import tree2

import sys, os, getopt

#sys.setrecursionlimit(20000)

#deventshort = {'S':'S', 'D':'D', 'H':'T'}
orinhxcommsep=':'
cleannhxcommsep='/'

def parse_branch_annot(annot):
	levt = []
	for elt in annot.split(cleannhxcommsep):
		if not elt=='&&NHX':
			evtype, evloc = elt.split('=')
			if evloc!='N':
				if (elt[0] in ('S', 'D')):
					levt.append((evtype, evloc, ''))
				elif elt[0]=='H':
					don, rec = evloc.split('@')[1:]
					levt.append(('T', rec, don))
	return levt
	


def parseMultiNHXRec(nfinnhx):
	lgt = []
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
