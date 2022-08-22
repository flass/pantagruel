#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys, os
import glob

nfin = sys.argv[1]
famprefix = sys.argv[2]
dirout = sys.argv[3]
padlen = int(sys.argv[4])
if len(sys.argv)>5: 
	writeseq = bool(int(sys.argv[5]))
else:
	writeseq = True
if len(sys.argv)>6: 
	discardsingle = bool(int(sys.argv[6]))
else:
	discardsingle = False
	
dseq2prevfam = {}
dprevfam2reprseq = {}
maxprevfamn = 0
if len(sys.argv)>7: 
	nfprevdbidentseq = sys.argv[7]
	print "will pick family numbering based on previous clustering (more intensive task)"
	with open(nfprevdbidentseq, 'r') as fprevdbidentseq:
		for line in fprevdbidentseq:
			fam, prot = line.rstrip('\n').split('\t')
			dseq2prevfam[prot] = fam
			if not fam in dprevfam2reprseq: dprevfam2reprseq[fam] = prot
			famn = int(fam.split(famprefix)[1]) # assumes the same prefix has been used! should throw an error if not
			if famn > maxprevfamn: maxprevfamn = famn # a bit over the top as the fam ids should come sequentially, but this is more robust

lvar = []
for var in ['nfin', 'famprefix', 'dirout', 'padlen', 'writeseq', 'discardsingle', 'nfprevdbidentseq']:
	lvar.append("%s = %s"%(var, repr(eval(var))))
print ' ; '.join(lvar)

if writeseq:
	if not os.path.exists(dirout):
		os.mkdir(dirout)
	else:
		raise IOError, "ouput directory '%s' already exists"%dirout
	

def idfam(nfam=-1, famprefix=famprefix, padlen=padlen):
	nfam += 1
	return (famprefix+str(nfam).zfill(padlen), nfam)
	
idfam0, nfam0 = idfam()
if writeseq: fout0 = open("%s/%s.fasta"%(dirout, idfam0), 'w')	# family with id # PREFIX000000 is reserved for ORFan sequences
ftabout = open("%s.tab"%dirout, 'w')

def incrementfam(nfam, lseqinfam, seqbuffer, prevfam=None, dprevfam2reprseq=None):
	if len(lseqinfam)>1:
		if prevfam:
			idfamn = prevfam
			if dprevfam2reprseq:
				reprseq = dprevfam2reprseq[prevfam]
				irepr = lseqinfam.index(reprseq)
				lseqinfam = [lseqinfam[irepr]]+lseqinfam[:irepr]+lseqinfam[(irepr  + 1):]
		else:
			# generate new id
			idfamn, nfam = idfam(nfam)
		if writeseq:
			# open, fill and close family file
			with open("%s/%s.fasta"%(dirout, idfamn), 'w') as fout: fout.write(seqbuffer)
	else:
		idfamn = idfam0
		if writeseq:
			fout0.write(seqbuffer)
	if not (idfamn==idfam0 and discardsingle):
		ftabout.writelines(["%s\t%s\n"%(idfamn, seqname) for seqname in lseqinfam])
	# update/reset
	lseqinfam = []
	seqbuffer = ''
	return (nfam, lseqinfam, seqbuffer)

nfam = maxprevfamn
lseqinfam = []
seqbuffer = ''

lnfinn = glob.glob(nfin) + glob.glob(nfin+'.[0-9]*') # left glob search looks for output from MMseqs v7 and prior version; right glob search looks for output from  MMseqs v8 and later
for nfinn in lnfinn:
	with open(nfinn, 'r') as finn:
		for line in finn:
			if line.startswith('\0'):
				# write sequence buffer for n-1 family and reset it
				if dseq2prevfam and len(lseqinfam)>1:
					# search previous clustering for an already assigned family
					# (do not search for singletons, bound to be all in the same 'family' PREFIX000000)
					prevfams = set([])
					for seqname in lseqinfam:
						prevfam = dseq2prevfam.get(seqname)
						if prevfam: prevfams.add(prevfam)
					if prevfams:
						if len(prevfams)>1:
							raise IndexError, "too many family ids from preious clustering associated to proteins that are members of a single family in new clustering:\nprevious clustering families: "+' '.join(prevfams)+'\n'+'\n'.join(lseqinfam[:2])+'\n...\n'+'\n'.join(lseqinfam[-2:])
						else:
							nfam, lseqinfam, seqbuffer = incrementfam(nfam, lseqinfam, seqbuffer, prevfam=list(prevfams)[0], dprevfam2reprseq=dprevfam2reprseq)
				else:
					nfam, lseqinfam, seqbuffer = incrementfam(nfam, lseqinfam, seqbuffer)
				l = line.lstrip('\0')
			else:
				l = line
			if l.startswith('>'):
				seqname = l.split(' ', 1)[0].lstrip('>')
				lseqinfam.append(seqname)
			seqbuffer += l
		
