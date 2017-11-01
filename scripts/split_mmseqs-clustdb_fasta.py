#!/usr/bin/python
import sys, os
nfin = sys.argv[1]
famprefix = sys.argv[2]
if len(sys.argv)>3: padlen = sys.argv[3]
else: padlen = 6

dirout = "%s_fasta"%nfin
if not os.path.exists(dirout):
	os.mkdir(dirout)
else:
	raise IOError, "ouput directory '%s' already exists"%dirout
	

def idfam(nfam):
	return famprefix+str(nfam).zfill(padlen)
	
idfam0 = idfam(0)
fout0 = open("%s/%s.fasta"%(dirout, idfam0), 'w')	# family with id # PREFIX000000 is reserved for ORFan sequences
ftabout = open("%s.tab"%dirout, 'w')

def incrementfam(nfam, lseqinfam, seqbuffer):
	if len(lseqinfam)>1:
		# generate new id
		idfamn = idfam(nfam)
		nfam += 1
		# open, fill and close family file
		with open("%s/%s.fasta"%(dirout, idfamn), 'w') as fout: fout.write(seqbuffer)
		ftabout
	else:
		idfamn = idfam0
		fout0.write(seqbuffer)
	ftabout.writelines(["%s\t%s\n"%(idfamn, seqname) for seqname in lseqinfam])
	# update/reset
	lseqinfam = []
	seqbuffer = ''
	return (nfam, lseqinfam, seqbuffer)

nfam = 1
lseqinfam = []
seqbuffer = ''

with open(nfin, 'r') as fin:
	for line in fin:
		if line.startswith('\0'):
			# write buffer for previous family and reset it
			nfam, lseqinfam, seqbuffer = incrementfam(nfam, lseqinfam, seqbuffer)
			l = line.lstrip('\0')
		else:
			l = line
		if l.startswith('>'):
			seqname = l.split()[0].lstrip('>')
			lseqinfam.append(seqname)
		seqbuffer += l
		
