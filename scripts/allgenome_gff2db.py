#!/usr/bin/python

import sys
import os
import gzip
import re

nfldirassemb = sys.argv[1]
dirout = sys.argv[2]
nftaxnamesdump = sys.argv[3]

if os.path.exists(dirout):
	if os.path.isfile(dirout):
		raise ValueError, "Cannot create folder '%s' as it already exist as a file"%(dirout)
else:
	os.mkdir(dirout)

with open(nfldirassemb, 'r') as fldirassemb:
	ldirassemb = [line.rstrip('\n') for line in fldirassemb.readlines()]

# link to NCBI Taxonomy database
dtaxid2sciname = {}
with open(nftaxnamesdump, 'r') as ftaxnamesdump:
	for line in ftaxnamesdump:
		lsp = line.rstrip('\n').split('\t')
		dtaxid2sciname[int(lsp[0])] = lsp[1]

daliasrepli = {'ANONYMOUS':'chromosome'}

assembpat = re.compile('(GCF_[^\._]+\.[0-9])_(.+)')

# prepare annoatation db output files
dfoutheaders = { \
'proteins':['protein_id', 'genomic_accession', 'gene_begin', 'gene_end', 'gene_strand', 'gene_product', 'genbank_cds_id'], \
'ribosomes':['rrna_id', 'genomic_accession', 'gene_begin', 'gene_end', 'gene_strand', 'gene_product'], \
'replicons':['assembly_accession', 'assembly_name', 'genomic_accession', 'replicon_name', 'replicon_type', 'size', 'tax_id', 'strain'] \
}
if dtaxid2sciname: dfoutheaders['replicons'].append('taxon_name')

dfout = {}
fouttags = ['proteins', 'ribosomes', 'replicons']
for fouttag in fouttags:
	dfout[fouttag] = open(os.path.join(dirout, "all%s_info.tab"%(fouttag) ), 'w')
	dfout[fouttag].write('\t'.join(dfoutheaders[fouttag])+'\n')

for dirassemb in ldirassemb:
	# extract assembly acc and name
	assembsearch = assembpat.search(os.path.basename(dirassemb))
	assacc, assname = assembsearch.groups()
	descprevgene = {}
	# parse CDS fasta file
	fcds = gzip.open("%s/%s_cds_from_genomic.fna.gz"%(dirassemb, os.path.basename(dirassemb)), 'rb')
	# register unambiguously the exact naming of the CDS sequence in the corresponding extracted CDS sequence file
	dgenbankcdsids = {}
	for line in fcds:
		if line.startswith('>'):
			genbankcdsid, cdstags = line.strip('>\n').split(' ', 1)
			dcdstags = dict(cdstag.split('=', 1) for cdstag in cdstags.strip('[]').split('] ['))
			dgenbankcdsids[(dcdstags['locus_tag'], dcdstags.get('protein_id', ''))] = genbankcdsid
	fcds.close()
	# parse GFF file
	fgff = gzip.open("%s/%s_genomic.gff.gz"%(dirassemb, os.path.basename(dirassemb)), 'rb')
	
	print "parse '%s'"%dirassemb
	currseqreg = ''
	cdslineout = None
	#~ numcds = 0
	nprintmullicds = 0
	for line in fgff:
		if not line.startswith('#'):
			lsp = line.rstrip('\n').split('\t')
			seqreg = lsp[0]
			annottype = lsp[2]
			beg = lsp[3]
			end = lsp[4]
			strand = lsp[6]
			desc = dict(d.split('=', 1) for d in lsp[8].split(';'))
			if lsp[2]=='region':
				if seqreg!=currseqreg and beg=='1':
					# new replicon
					currseqreg = seqreg
					dbxref = dict(d.split(':') for d in desc['Dbxref'].split(','))
					taxid = dbxref['taxon']
					strain = desc.get('strain', '')
					repliname = desc.get('Name', '')
					replitype = desc.get('genome', '')
					replineout = [assacc, assname, seqreg, daliasrepli.get(repliname, repliname), replitype, end, taxid, strain]
					if dtaxid2sciname: replineout.append(dtaxid2sciname[int(taxid)])
					dfout['replicons'].write('\t'.join(replineout)+'\n')
			elif lsp[2]=='gene':
				descprevgene = desc
				# need to buffer extraction of CDS entries as several lines can relate to the same gene in case of CDS in several segments (due  to introns, framshifts, ...)
				ngenecds = 0
				buffbeg = []
				buffend = []
				lastprodid = ''
			elif lsp[2] in ['CDS', 'rRNA']:
				if lsp[2]=="CDS":
					if ngenecds==0 and cdslineout:
						# write out last CDS info (potential synthesis of  several lines)
						dfout['proteins'].write('\t'.join(cdslineout)+'\n')
						cdslineout = None
					ngenecds += 1
					buffbeg.append(beg)
					buffend.append(end)
					productid = desc.get('protein_id', '')
					locustag = descprevgene['locus_tag']
					if ngenecds > 1:
						if nprintmullicds==0: print "multiline CDS feature:",
						print productid,
						nprintmullicds +=1
						if lastprodid!=productid: raise IndexError, "multiline CDS feature not pointing at the same protein: %s and %s"%(lastprodid, productid)
					cdslineout = [productid, seqreg, locustag, min(buffbeg), max(buffend), strand, desc.get('product', '')]
					# generate unique CDS id
					#~ numcds = int(desc['ID'].lstrip('cds'))+1 #feature numbering is 0-based; we want a 1-based numbering in output.
					#~ if productid: genbankcdsid = "lcl|%s_cds_%s_%d"%(seqreg, productid, numcds)
					#~ else: genbankcdsid = "lcl|%s_cds_%d"%(seqreg, numcds)
					genbankcdsid = dgenbankcdsids[(locustag, productid)]
					cdslineout.append(genbankcdsid)
				if lsp[2]=="rRNA":
					productid = desc.get('ID', '')
					riblineout = [productid, seqreg, locustag, beg, end, strand, desc.get('product', '')]
					dfout['ribosomes'].write('\t'.join(riblineout)+'\n')
				lastprodid = productid
	else:
		if cdslineout:
			# write out last CDS info (potential synthesis of  several lines); flush buffer
			dfout['proteins'].write('\t'.join(cdslineout)+'\n')
	if nprintmullicds>0: print	'\n'

for fouttag in fouttags:
	dfout[fouttag].close()
