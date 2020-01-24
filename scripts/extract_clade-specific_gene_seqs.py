#! /usr/bin/pyton2.7
import os, sys, glob, gzip
from Bio import SeqRecord, SeqIO
import sqlite3

joinlocusdist = 2000

#dirbgass = '/pantagruel_databases/choleraclub315plus7/00.input_data/genbank-format_assemblies'
#dout = '/pantagruel_databases/choleraclub315plus7/clade-specific_genefams/clade_specific_sequences'

nfcladespegenetab = sys.argv[1]
dirbgass = sys.argv[2]
dout = sys.argv[3]
nfdb = sys.argv[4]

dbcon = sqlite3.connect(nfdb)
dbcur = dbcon.cursor() 

dcode2refass = {}
code2assq = "select code, assembly_id, assembly_name from assemblies;"
dbcur.execute(code2assq)
lcoderefass = dbcur.fetchall()
for coderefass in lcoderefass:
	dcode2refass[coderefass[0]] = '_'.join(coderefass[1:])

dbcon.close()
	
#dclade2refass = {'clade236': ('VIBCHO30', '48611_A01.contigs_pacbio.1_Vibrio_cholerae_SMIC_67_01_ext2'),
#				 'clade152': ('VIBCHO400', 'SAMEA4384428.contigs.1_Vibrio_cholerae_noname400'),
#				}
##assemb = '48611_A01.contigs_pacbio.1_Vibrio_cholerae_SMIC_67_01_ext2'
##code = 'VIBCHO30'
##cla = 'clade236'
##claspeset = 'relax0.0_vscontrastclade_specific_pres_genes'
##specoords = [(1456822,1493162), (3206243, 3250862), (3313853, 3314119), (3367113, 3371021)]
#dclaspeset2coords = {'clade236': {'relax0.0_vscontrastclade_specific_pres_genes':[(0,1456822,1493162), (0,3206243, 3250862), (0,3313853, 3314119), (0,3367113, 3371021)],
#                                  'relax0.05_vscontrastclade_specific_pres_genes':[(0,1125281,1137921), (0,1456822,1493162), (0,1540035,1543530), (0,3206243, 3250862), (0,3313853, 3314119), (0,3367113, 3371021)]},
#					 'clade152': {'relax0.05_vscontrastclade_specific_pres_genes':[(1,86432,176408), (12,9957,18160), (15,63236,107467), ]}
#					}

claspeset = os.path.basename(nfcladespegenetab).rsplit('.', 1)[0]
fcladespegenetab = open(nfcladespegenetab, 'r')
dclaspeset2coords = {}
clablockcount = 0
spegeneclusters = []
header = None
cla = None
currspegenecluster = []

for line in fcladespegenetab:
	lsp = line.strip('\n# ').split('\t')
#	print lsp, clablockcount
	if line.startswith('#'):
		if clablockcount==0:
			if cla:
				if currspegenecluster:
					spegeneclusters.append(tuple(currspegenecluster))
					currspegenecluster = []
#				print cla, clablockcount, spegeneclusters, lsp
				dclaspeset2coords.setdefault(cla, {})[claspeset] = spegeneclusters
				spegeneclusters = []
			try:
				cla, cladename = lsp[0].split()
			except ValueError:
				continue
			if not cla.startswith('clade'):
				cla = None
				continue
#		elif clablockcount==1:
#			clacodes = lsp[7].split(',')
#		elif clablockcount==2:
#			bgclacodes = lsp[7].split(',')
		elif clablockcount==3:
			# no specific gene found for this clade
			clablockcount = 0
			currspegenecluster = []
			continue
		clablockcount += 1
	elif clablockcount>0:
		header = lsp
		clablockcount = 0
		ibeg = header.index('cds_begin')
		iend = header.index('cds_end')
		icontig = header.index('genomic_accession')
	else:
		gbeg = int(lsp[ibeg])-1
		gend = int(lsp[iend])-1
		contigid = lsp[icontig]
		if not currspegenecluster:
			currspegenecluster = [contigid, gbeg, gend]
		else:
			if contigid==currspegenecluster[0] and gbeg < currspegenecluster[2]+joinlocusdist:
				currspegenecluster[2] = gend
			else:
				spegeneclusters.append(tuple(currspegenecluster))
				currspegenecluster = [contigid, gbeg, gend]

# store last values
if cla:
	if currspegenecluster: spegeneclusters.append(tuple(currspegenecluster))
	dclaspeset2coords.setdefault(cla, {})[claspeset] = spegeneclusters

fcladespegenetab.close()
#print dclaspeset2coords



for cla in dclaspeset2coords:
	print cla
	for claspeset, specoords in dclaspeset2coords[cla].iteritems():
		if specoords:
			coderepr = specoords[0][0].rsplit('_', 1)[0]
			assemb = dcode2refass[coderepr]
		else:
			continue
		with gzip.open(os.path.join(dirbgass, assemb, '%s_genomic.gbff.gz'%assemb), 'rb') as fgbff:
			dseqrecid2seqrec = {seqrec.id:seqrec for seqrec in SeqIO.parse(fgbff, format='genbank')}
		nfout = os.path.join(dout, '%s-%s-%s-%s.fasta'%(cla, claspeset, coderepr, assemb))
		sperecs = []
		for i, specoord in enumerate(specoords):
#			print specoord
			contigid, locusbeg, locusend = specoord
			sperec = dseqrecid2seqrec[contigid][(locusbeg):(locusend)]
			sperec.id = sperec.id+'-spegenecluster%d-%d-%d'%(i, locusbeg, locusend)
			sperecs.append(sperec)
		SeqIO.write(sperecs, handle=nfout, format='fasta')

