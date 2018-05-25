#!/usr/bin/python

### general purpose functions

def mean(seq, ignoreNull=True):
	l = [float(k) for k in seq if ((k is not None) or (not ignoreNull))]
	if not l: return None
	return sum(l)/len(l)

def median(seq, ignoreNull=True):
	l = [k for k in seq if ((k is not None) or (not ignoreNull))]
	l.sort()
	L = len(l)
	if not l: return None
	nmed = L/2
	if L%2==1: return float(l[nmed])
	else: return (float(l[nmed-1])+float(l[nmed]))/2

def var(seq, correct=1):
	m = mean(seq)
	if m is None: return None
	n = len([k for k in seq if (k is not None)])-correct
	if n < 1:
		return 'inf'		
	else:
		d = [(float(x) - m)**2 for x in seq]
		return sum(d)/n

def findSeqRecordIndexesFromSeqNames(aln, seqnames):
	if isinstance(seqnames, list): return [k for k,seq in enumerate(aln) if seq.id in seqnames]
	else:
		for k,seq in enumerate(aln):
			if seq.id in seqnames: return k
