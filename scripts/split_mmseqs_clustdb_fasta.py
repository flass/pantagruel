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
	
lvar = []
for var in ['nfin', 'famprefix', 'dirout', 'padlen', 'writeseq', 'discardsingle']:
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

def incrementfam(nfam, lseqinfam, seqbuffer):
	if len(lseqinfam)>1:
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

nfam = 0
lseqinfam = []
seqbuffer = ''

lnfinn = glob.glob(nfin) + glob.glob(nfin+'.[0-9]*') # left glob search looks for output from MMseqs v7 and prior version; right glob search lloks for output from  MMseqs v8 and later
for nfinn in lnfinn:
	with open(nfinn, 'r') as finn:
		for line in finn:
			if line.startswith('\0'):
				# write buffer for previous family and reset it
				nfam, lseqinfam, seqbuffer = incrementfam(nfam, lseqinfam, seqbuffer)
				l = line.lstrip('\0')
			else:
				l = line
			if l.startswith('>'):
				seqname = l.split(' ', 1)[0].lstrip('>')
				lseqinfam.append(seqname)
			seqbuffer += l
		
