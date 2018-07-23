#!/usr/bin/python
import os, sys
import sqlite3

def writefamrowout(currphyloprof, curfam, lgenomecodes, dfoutmatortho, lcurrt):
	if curfam[1] is not None: scurfam = '%s-%d'%curfam
	else: scurfam = curfam[0]
	k = sum([int(bool(n)) for n in currphyloprof.values()])
	if k==0: raise ValueError, "no gene assigned in fam %s\nprevious rows were:\n%s"%(scurfam, repr(lcurrt))
	elif k==1: foutmatortho = dfoutmatortho['singletons']
	else: foutmatortho = dfoutmatortho['no-singletons']
	foutmatortho.write('\t'.join([scurfam]+[str(currphyloprof[gcode]) for gcode in lgenomecodes])+'\n')
	lcurrt = []
	return (k==1)

if len(sys.argv) < 2:
	print "Missing arguments!\nUsage:\n%s /path/to/pantagruel/sqliteDBfile /path/to/output/[file_prefix] [/path/to/restricted_genome_code_list] [ortholog_collection_id]"%os.path.basename(sys.argv[0])
	sys.exit(2)

nfsqldb = sys.argv[1]
nfoutrad = sys.argv[2].rstrip('/')
if os.path.dirname(nfoutrad)==nfoutrad:
	nfoutrad = os.path.join((nfoutrad, 'pantagruel_orthologs'))

if len(sys.argv) > 3:
	nffocusgenomecodes = sys.argv[3]
else:
	nffocusgenomecodes = None
	

if len(sys.argv) > 4:
	ortcolid = int(sys.argv[3])
else:
	ortcolid = None


dbcon = sqlite3.connect(nfsqldb)
dbcur = dbcon.cursor()

if nffocusgenomecodes:
	with open(nffocusgenomecodes, 'r') as ffocusgenomecodes:
		lfocusgenomecodes = [(line.strip(' \n'),) for line in ffocusgenomecodes]
	dbcur.execute("CREATE TEMP TABLE focus_genome_codes (code VARCHAR(10) NOT NULL);")
	dbcur.executemany("INSERT INTO focus_genome_codes VALUES (?);", lfocusgenomecodes)
	IJfocus = "INNER JOIN focus_genome_codes USING (code)"
	lgenomecodes = sorted([tcode[0] for tcode in lfocusgenomecodes])
else:
	dbcur.execute("SELECT code FROM assemblies ORDER BY code;")
	lgenomecodes = [tcode[0] for tcode in dbcur]
	IJfocus = ""

print "building matrix of gene presence / absence for %d genomes"%len(lgenomecodes)


dbcur.execute("""
CREATE TEMP TABLE cds_fam_code AS 
SELECT cds_code, gene_family_id, code
FROM coding_sequences 
INNER JOIN replicons USING (genomic_accession) 
INNER JOIN assemblies USING (assembly_id) %s 
WHERE gene_family_id IS NOT NULL AND gene_family_id!='RHIZOC000000'
ORDER BY gene_family_id, code;"""%(IJfocus))

dbcur.execute("SELECT count(*) FROM cds_fam_code;")
print "examining a total of %d CDSs with non-ORFan family assignment"%(dbcur.fetchone()[0])

dfoutmatortho = {}
for s in ['singletons', 'no-singletons']:
	nfoutmatortho = nfoutrad+"_genome_counts.%s.mat"%s
	foutmatortho = open(nfoutmatortho, 'w')
	dfoutmatortho[s]  = foutmatortho
	foutmatortho.write('\t'.join(['']+lgenomecodes)+'\n')

colWC = "WHERE ortholog_col_id=%d"%ortcolid if not (ortcolid is None) else ""
dbcur.execute("""
SELECT gene_family_id, og_id, code, count(cds_code)
FROM cds_fam_code
LEFT JOIN orthologous_groups using (cds_code, gene_family_id) %s
GROUP BY gene_family_id, og_id, code
ORDER BY gene_family_id, og_id;"""%colWC)

lnoorthofams = []
lorthofams = []
northofamogs = 0
northosing = 0
northonosing = 0
curfam = None
curfamog = None
lcurrt = []
currphyloprof = None
for tfamcodeogn in dbcur:
	lcurrt.append(tfamcodeogn)
	fam, ogid, code, n = tfamcodeogn
	if fam != curfam:
		curfam = fam
		if ogid is not None: lorthofams.append(fam)
		else:
			lnoorthofams.append(fam)
	if (fam, ogid) != curfamog:
		if curfamog:
			singleton = writefamrowout(currphyloprof, curfamog, lgenomecodes, dfoutmatortho, lcurrt)
			if singleton: northosing += 1
			else: northonosing += 1
		curfamog = (fam, ogid)
		currphyloprof = {gcode:0 for gcode in lgenomecodes}
		if ogid: northofamogs += 1
	currphyloprof[str(code)] = n
else:
	singleton = writefamrowout(currphyloprof, curfamog, lgenomecodes, dfoutmatortho, lcurrt)
	if singleton: northosing += 1
	else: northonosing += 1

nfoutlnoortho = nfoutrad+"_families_no-orthologs"
with open(nfoutlnoortho, 'w') as foutlnoortho:
	foutlnoortho.write('\n'.join(lnoorthofams)+'\n')
print "%d families not covererd by orthology classification (means no evolution scenario was inferred for these families)"%len(lnoorthofams)

nfoutlortho = nfoutrad+"_families_with_orthologs"
with open(nfoutlortho, 'w') as foutlortho:
	foutlortho.write('\n'.join(lorthofams)+'\n')
print "%d families covererd by orthology classification into a total of %d orthologous groups"%(len(lorthofams), northofamogs)

print "these toalize %d families with unique representative in the dataset (singletons) and %d others [total: %d]"%(northosing, northonosing, northosing+northonosing)

dbcon.close()
for s in ['singletons', 'no-singletons']:
	dfoutmatortho[s].close()
