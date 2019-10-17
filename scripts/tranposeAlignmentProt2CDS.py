#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""transpose a CDS into the positions of the corresponding aligned protein;
   assumes no indels, only mismatches and possibly shortenned sequences.
   THIS PROGRAM DOES NOT VERIFY THAT THE CDS TRANSLATES INTO THE PROTEIN.
   this is used as a fallback when pal2nal.pl fails, for instance due to 
   use of non-standard amino-acids like selenocystein, or the premature 
   termination of one of the sequences (typicaly the lack of the terminal 
   segment of the CDS not covered by genomic DNA record, while full(er) 
   protein sequence was charaterized)
   NB: late start of the CDS relative to the protein will lead to abberant alignment!!!
"""
from Bio import Alphabet, Seq, SeqIO
import sys, os

gap = '-'
ab = Alphabet.generic_dna
gAB = Alphabet.Gapped(Alphabet.generic_protein, gap)

nfcdsseq = os.path.abspath( sys.argv[1] )
nfprtali = os.path.abspath( sys.argv[2] )
nfout    = os.path.abspath( sys.argv[3] )

itcdsseq = SeqIO.parse(nfcdsseq, format='fasta', alphabet=ab)
itprtali = SeqIO.parse(nfprtali, format='fasta', alphabet=gAB)
fout = open(nfout, 'w')

for cdsseqrec in itcdsseq:
	prtalirec = itprtali.next()
	fout.write('>%s\n'%prtalirec.description)
	i = 0
	j = 0
	strcds = str(cdsseqrec.seq)
	strprot = str(prtalirec.seq)
	strout = ''
	C = len(strcds)
	P = len(strprot)
	while ((i+1)*3 <= C) and (j < P):
		codon = strcds[(i*3):((i+1)*3)]
		aa = strprot[j]
		i += 1 ; j += 1
		while aa==gap and (j < P):
			strout += gap*3
			aa = strprot[j]
			j += 1
		if j < P:
			# did not reach the end of line in the protein alignment; can add the sequence codon
			strout += codon
	while j < P-1:
		# did not reach the end of line in the protein alignment; needs padding with gaps for a clean result
		strout += gap*3 ; j += 1
	fout.write(strout+'\n')

fout.close()
print "CDS alignment generated at: '%s'"%(nfout)
print "WARNING: --- correctness of the alignment not guaranted, please verify ---"
