#!/usr/bin/python
import sys, os
nfin = sys.argv[1]
famprefix = sys.argv[2]
dirout = sys.argv[3]
padlen = int(sys.argv[4])
if len(sys.argv)>5: 
	writeseq = True
else:
	writeseq = bool(int(sys.argv[5]))
if len(sys.argv)>6: 
	discardsingle = False
else:
	discardsingle = bool(int(sys.argv[6]))
	
for var in ['nfin', 'famprefix', 'dirout', 'padlen', 'writeseq', 'discardsingle']:
	print "%s = %s"%(var, repr(eval(var))),
print ''

if writeseq and (not os.path.exists(dirout)):
	os.mkdir(dirout)
else:
	raise IOError, "ouput directory '%s' already exists"%dirout
	

def idfam(nfam, famprefix=famprefix, padlen=padlen):
	return famprefix+str(nfam).zfill(padlen)
	
idfam0 = idfam(0)
if writeseq: fout0 = open("%s/%s.fasta"%(dirout, idfam0), 'w')	# family with id # PREFIX000000 is reserved for ORFan sequences
ftabout = open("%s.tab"%dirout, 'w')

def incrementfam(nfam, lseqinfam, seqbuffer):
	if len(lseqinfam)>1:
		# generate new id
		idfamn = idfam(nfam)
		nfam += 1
		if writeseq:
			# open, fill and close family file
			with open("%s/%s.fasta"%(dirout, idfamn), 'w') as fout: fout.write(seqbuffer)
	else:
		idfamn = idfam0
		if writeseq:
			fout0.write(seqbuffer)
	if not (idfam==idfam0 and discardsingle):
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
		
