#!/usr/bin/python

import sys
import os
import gzip
import re

def parseGFFline(line):
	#~ print line
	lsp = line.rstrip('\n').split('\t')
	seqreg = lsp[0]
	annottype = lsp[2]
	beg = lsp[3]
	end = lsp[4]
	strand = lsp[6]
	desc = dict(d.split('=', 1) for d in lsp[8].split(';'))
	return (seqreg, annottype, beg, end, strand, desc)


nfldirassemb = sys.argv[1]
dirout = sys.argv[2]
dirncbitax = sys.argv[3]

if os.path.exists(dirout):
	if os.path.isfile(dirout):
		raise ValueError, "Cannot create folder '%s' as it already exist as a file"%(dirout)
else:
	os.mkdir(dirout)

with open(nfldirassemb, 'r') as fldirassemb:
	ldirassemb = [line.rstrip('\n') for line in fldirassemb.readlines()]

# link to NCBI Taxonomy database
dmergedtaxid = {}
with open("%s/merged.dmp"%dirncbitax, 'r') as ftaxmergedump:
	for line in ftaxmergedump:
		lsp = [f.strip('\t') for f in line.rstrip('|\n').split('|')]
		dmergedtaxid[int(lsp[0])] = int(lsp[1])

dtaxid2sciname = {}
with open("%s/names.dmp"%dirncbitax, 'r') as ftaxnamesdump:
	for line in ftaxnamesdump:
		lsp = [f.strip('\t') for f in line.rstrip('|\n').split('|')]
		if lsp[-1]=='scientific name':
			dtaxid2sciname[int(lsp[0])] = lsp[1]

	
daliasrepli = {'ANONYMOUS':'chromosome'}

assembpatgenbank = re.compile('(GC[AF]_[^\._]+\.[0-9])_(.+)')
assembpat = re.compile('(.+?\.[0-9])_(.+)')


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
	# first scans the file for region and (pseudo)gene features
	dgeneranges = {}
	dgenedesc = {}
	#~ dregions = {}
	for line in fgff:
		if line.startswith('#'): continue
		elif line.startswith('>'): break
		#~ print line
		seqreg, annottype, beg, end, strand, desc = parseGFFline(line)
		if annottype=='region':
			if seqreg!=currseqreg and beg=='1':
				# new replicon
				#~ dregions[desc['ID']] = desc
				currseqreg = seqreg.split('|')[-1]	# to account for Prokka-style annotation that prepends 'gnl|Sequencing_Centre' to the region_id
				dbxref = dict(d.split(':') for d in desc['Dbxref'].split(','))
				taxid = dbxref['taxon']
				strain = desc.get('strain', '')
				repliname = desc.get('Name', '')
				replitype = desc.get('genome', '')
				replineout = [assacc, assname, seqreg, daliasrepli.get(repliname, repliname), replitype, end, taxid, strain]
				tid = int(taxid)
				if dtaxid2sciname: replineout.append(dtaxid2sciname[dmergedtaxid.get(tid, tid)])
				dfout['replicons'].write('\t'.join(replineout)+'\n')
		elif annottype in ['gene', 'pseudogene']:
			#~ descprevgene = desc
			dgeneranges[desc['ID']] = []
			dgenedesc[desc['ID']] = desc
	
	# need to buffer extraction of CDS entries as several lines can relate to the same gene in case of CDS in several segments (due  to introns, framshifts, ...)
	buffbeg = []
	buffend = []
	lastprodid = ''
	# resume reading the file from start
	fgff.seek(0)
	for line in fgff:
		if line.startswith('#'): continue
		elif line.startswith('>'): break
		#~ print line
		seqreg, annottype, beg, end, strand, desc = parseGFFline(line)
		if annottype in ['CDS', 'rRNA']:
			if annottype=="CDS":
				parentgene = desc['Parent']
				if len(dgeneranges[parentgene])==0 and cdslineout:
					# write out last CDS info (potential synthesis of  several lines)
					dfout['proteins'].write('\t'.join(cdslineout)+'\n')
					cdslineout = None
					lastprodid = ''
				productid = desc.get('protein_id', '').split('|')[-1]	# to account for Prokka-style annotation that prepends 'gnl|Sequencing_Centre' to the protein_id
				dgeneranges[parentgene].append((beg, end, productid))
				begs, ends, products = zip(*dgeneranges[parentgene])
				locustag = dgenedesc[parentgene]['locus_tag']
				if len(dgeneranges[parentgene]) > 1:
					if nprintmullicds==0: print "multiline CDS feature:",
					print productid,
					nprintmullicds +=1
					if products[0]!=productid: raise IndexError, "multiline CDS feature not pointing at the same protein: %s and %s"%(products[0], productid)
				cdslineout = [productid, seqreg, locustag, min(begs), max(ends), strand, desc.get('product', '')]
				# generate unique CDS id
				#~ numcds = int(desc['ID'].lstrip('cds'))+1 #feature numbering is 0-based; we want a 1-based numbering in output.
				#~ if productid: genbankcdsid = "lcl|%s_cds_%s_%d"%(seqreg, productid, numcds)
				#~ else: genbankcdsid = "lcl|%s_cds_%d"%(seqreg, numcds)
				genbankcdsid = dgenbankcdsids[(locustag, productid)]
				cdslineout.append(genbankcdsid)
			if annottype=="rRNA":
				productid = desc.get('ID', '')
				riblineout = [productid, seqreg, locustag, beg, end, strand, desc.get('product', '')]
				dfout['ribosomes'].write('\t'.join(riblineout)+'\n')
	if cdslineout:
		# write out last CDS info (potential synthesis of  several lines); flush buffer
		dfout['proteins'].write('\t'.join(cdslineout)+'\n')
	if nprintmullicds>0: print	'\n'

for fouttag in fouttags:
	dfout[fouttag].close()
