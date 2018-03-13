#!/usr/bin/python

import sys, os

minchrlen = 2000000

nfgffin = sys.argv[1]
nfgffout = sys.argv[2]
nfstraininfo = sys.argv[3]
nfrawassembseq = sys.argv[4]
assemblervers = sys.argv[5]

dstraininfo = {}
with open(nfstraininfo, 'r') as fstraininfo:
	header = fstraininfo.readline().rstrip('\n').split('\t')
	for line in fstraininfo:
		lsp = line.rstrip('\n').split('\t')
		loctagprefix = lsp[5]
		dstraininfo[loctagprefix] = dict(zip(header, lsp))

lcontignames = []
dcontigs = {}		
with open(nfrawassembseq, 'r') as frawassembseq:
	for line in frawassembseq:
		if line.startswith('>'):
			lsp = line.strip('>\n').split()
			if lsp[0].startswith('NODE'):
				# velvet/spades type
				lsp = lsp[0].split('_')
				cname = '_'.join(lsp[:2])
				dcontigs[cname] = {lsp[i]:lsp[i+1] for i  in range(1, len(lsp)/2)}					
			else:
				try:
					# Unicylcer type, name is just a numeral
					cname = 'contig'+str(int(lsp[0]))
				except ValueError:
					cname = lsp[0]
				dcontigs[cname] = dict(lsp[i].split('=') for i  in range(1, len(lsp)))
			lcontignames.append(cname)

fgffin = open(nfgffin, 'r')
fgffout = open(nfgffout, 'w')
lregions = []
lregionlens = []
currseqreg = None
nreg = 0
fastatime = False

for line in fgffin:
	if line.startswith('##'):
		if line.startswith('##sequence-region'):
				lsp = line.rstrip('\n').split()
				lregions.append(lsp[1])
				lregionlens.append(lsp[2:4])
		if line.startswith('##FASTA'):
			fastatime = True
	else:
		seqreg = line.split('\t')[0]
		if not fastatime and (seqreg != currseqreg):
			# get strain info
			loctagprefix = seqreg.split('|')[-1].split('_')[0]
			straininfo = dstraininfo[loctagprefix]
			# add organism info line
			if nreg==0:
				fgffout.write('##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s\n'%straininfo['taxid'])
			currseqreg = seqreg
			#~ print nreg, seqreg
			# add region info line
			iscircular = dcontigs[lcontignames[nreg]].get('circular', 'false')
			replitype = 'chromosome' if (nreg==0 and lregionlens[nreg][1]>minchrlen) else ('plasmid' if iscircular=='true' else 'contig')
			extfeat = ['ID=id%d'%nreg, 'Dbxref=taxon:%s'%straininfo['taxid'], 'Is_circular=%s'%iscircular, 'gbkey=Src', 'genome=%s'%replitype, 'mol_type=genomic DNA', 'strain=%s'%straininfo['strain']]
			cols = [seqreg, assemblervers, 'region', lregionlens[nreg][0], lregionlens[nreg][1], '.', '+', '.', ';'.join(extfeat)]
			regline = '\t'.join(cols)+'\n'
			fgffout.write(regline)
			print regline.rstrip('\n')
			nreg += 1
	fgffout.write(line)

fgffin.close()
fgffout.close()
