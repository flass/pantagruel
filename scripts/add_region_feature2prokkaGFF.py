#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os

minchrlen = 2000000
verbose = False

gproject = sys.argv[1]
nfgffin = sys.argv[2]
nfgffout = sys.argv[3]
nfstraininfo = sys.argv[4]
nfrawassembseq = sys.argv[5]
assemblervers = sys.argv[6]
#nfgproject2assemb = sys.argv[6]
if len(sys.argv)>7:
	if sys.argv[7]=='verbose':
		verbose=True

dstraininfo = {}
with open(nfstraininfo, 'r') as fstraininfo:
	header = fstraininfo.readline().rstrip('\n').split('\t')
#	loctagprefixfield = header.index('locus_tag_prefix')
	gprojfield = header.index('assembly_id')
	for line in fstraininfo:
		lsp = line.rstrip('\n').split('\t')
		gproj = lsp[gprojfield]
		dstraininfo[gproj] = dict(zip(header, lsp))

straininfo = dstraininfo[gproject]
print "extracted strain information:"
print '\t'.join(header)
print '\t'.join(dstraininfo[gproject])
#dgproject2assemb = {}
#with open(nfgproject2assemb, 'r') as fgproject2assemb:
#	header = fgproject2assemb.readline().rstrip('\n').split('\t')
#	gsourcefield = header.index('genome_source')
#	for line in fgproject2assemb:
#		lsp = line.rstrip('\n').split('\t')
#		gsource = lsp[gsourcefield]
#		dgproject2assemb[gsource] = dict(zip(header, lsp))
	

lcontignames = []
dcontigs = {}
dcontiglenids = {}
clen = 0
cname = ''
c = 0
with open(nfrawassembseq, 'r') as frawassembseq:
	for line in frawassembseq:
		if line.startswith('>'):
			if clen>0 and cname:
				# record the contig length fo sorting them by decreasing length (they should already be)
				dcontiglenids[cname] = (clen, c)
				# also record the original position so to break ties by keping the original order
				c -= 1		# count into negatives to use as a reverse sorting key
				clen = 0
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
				try:
					dcontigs[cname] = dict(lsp[i].split('=') for i  in range(1, len(lsp)))
				except ValueError:
					dcontigs[cname] = {}
			lcontignames.append(cname)
		else:
			clen += len(line.strip(' \n'))
	else:
		if clen>0 and cname:
			dcontiglenids[cname] = (clen, c)
			c -= 1		# count into negatives to use as a reverse sorting key
			clen = 0

# sort on decreasing contig length (and on ex-aequo, on the original order)
lcontignames.sort(key=lambda x: dcontiglenids[x], reverse=True)
#print lcontignames
			
fgffin = open(nfgffin, 'r')
if verbose: print "open input GFF: '%s'"%nfgffin
fgffout = open(nfgffout, 'w')
if verbose: print "open output GFF: '%s'"%nfgffout
lregions = []
#~ lregionlens = []
dregionlens = {}
currseqreg = None
nregin = 0
nregout = 0
fastatime = False
dgffcontigname2rawcontigname = {}
#loctagprefix = None
#
for line in fgffin:
	if line.startswith('##sequence-region'):
		lsp = line.rstrip('\n').split()
		seqreg = lsp[1]
		lregions.append(seqreg)
		dregionlens[seqreg] = tuple(lsp[2:4])
#			dgffcontigname2rawcontigname[lsp[1]] = lcontignames[nregin]
#			# check this is the same lentgh as length-ordered contigs from the original contig file
#			print nregin, repr(lsp), dcontiglenids[lcontignames[nregin]]
#			if not int(lsp[3])==dcontiglenids[lcontignames[nregin]][0]:
#				raise IndexError, "unmatched contig order:\nconting from original file: %s\n#%d 'sequence-region' annotation from Prokka GFF: %s"%(dcontiglenids[lcontignames[nregin]], nregin, repr(lsp))
		nregin += 1
#		if not loctagprefix:
#			# get strain info - just need to do that once
#			loctagprefix = seqreg.split()[0].split('|')[-1].rsplit('_', 1)[0]
#			try:
#				straininfo = dstraininfo[loctagprefix]
#			except KeyError, e:
#				if loctagprefix == 'chr':
#					# special case of annotation of the chromosome contig with a different tag
#					# believed to be a legacy behaviour of prokka (<= v1.11), cf. https://github.com/flass/pantagruel/issues/25
#					# just wait for the next contig
#					loctagprefix = None
#				else:
#					raise KeyError, "missing information on genome with locus_tag prefix '%s'; please complete info in file '%s'"%(loctagprefix, nfstraininfo)
#
	elif not line.startswith('##gff-version'):
		fgffout.write(line)
		break
	fgffout.write(line)

#if not loctagprefix:
#	raise IndexError, """could not find out what is the locus_tag_prefix string from parsing the '##sequence-region' fields 
#	                     in the input GFF file '%s' (excluding those only annotated as 'chr') 
#						 so won't be able to match it to the strain reference file '%s'"""%(nfgffin, nfstraininfo)

# sort on decreasing contig length
lregions.sort(key=lambda x: dregionlens[x], reverse=True)
#print lregions
if verbose: print "len(lregions)", len(lregions), "len(lcontignames)", len(lcontignames)
dgffcontigname2rawcontigname = dict(zip(lregions, lcontignames))
#print dgffcontigname2rawcontigname

for line in fgffin:
	# resume input GFF file exploration
	if line.startswith('##'):
		if line.startswith('##FASTA'):
			fastatime = True
	else:
		seqreg = line.split('\t')[0]
		if not fastatime and (seqreg != currseqreg):
			# add organism info line
			if nregout==0:
				fgffout.write('##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s\n'%straininfo['taxid'])
			currseqreg = seqreg
			#~ print nregout, seqreg
			# add region info line
			#~ iscircular = dcontigs[lcontignames[nregout]].get('circular', 'false')
			rawcontigname = dgffcontigname2rawcontigname[seqreg]
			iscircular = dcontigs[rawcontigname].get('circular', 'false')
			#~ replitype = 'chromosome' if (nregout==0 and int(lregionlens[nregout][1])>minchrlen) else ('plasmid' if iscircular=='true' else 'contig')
			replitype = 'chromosome' if (nregout==0 and int(dregionlens[seqreg][1])>minchrlen) else ('plasmid' if iscircular=='true' else 'contig')
			extfeat = ['ID=id%d'%nregout, 'Dbxref=taxon:%s'%straininfo['taxid'], 'Is_circular=%s'%iscircular, 'gbkey=Src', 'genome=%s'%replitype, 'mol_type=genomic DNA', 'strain=%s'%straininfo['strain']]
			#~ cols = [seqreg, assemblervers, 'region', lregionlens[nregout][0], lregionlens[nregout][1], '.', '+', '.', ';'.join(extfeat)]
			cols = [seqreg, assemblervers, 'region', dregionlens[seqreg][0], dregionlens[seqreg][1], '.', '+', '.', ';'.join(extfeat)]
			regline = '\t'.join(cols)+'\n'
			fgffout.write(regline)
			print regline.rstrip('\n')
			nregout += 1
	fgffout.write(line)

fgffin.close()
fgffout.close()
