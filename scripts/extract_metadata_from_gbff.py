#!/usr/bin/python
import sys, os, getopt
import gzip
import re

## constants

# time format matching
dmonalias = dict(zip(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], [str(n) for n in range(1,13)]))
ldatepats = [ \
re.compile("([0-9]{2})-([A-Z][a-z]{2})-([0-9]{4})"), \
re.compile("([0-9]{4})-([A-Z][a-z]{2})-([0-9]{2})"), \
re.compile("([A-Z][a-z]{2})-([0-9]{4})"), \
re.compile("([0-9]{4})-([0-9]{2})"), \
re.compile("([0-9]{4})") ]

qualifgroups =	[["air", "water", "soil", "environment", "sink", "garden", "wild", "wastewater", "environmental"], \
				 ["food", "meat", "beef", "turkey", "butter", "meal", "broiler", "egg", "rice", "tuna", "steamed", "foodborne"], \
				 ["carcass", "abbatoir"], \
				 ["stool", "feces", "cecal", "cecum", "fecal", "duodenum", "perirectal", "gut", "rectum", "intestinal"], \
				 ["alfalfa", "apple", "lettuce", "tomatoes", "root", "marjoram", "vegetal"], \
				 ["endosymbiotic", "endosymbiont"], \
				 ["symbiotic", "symbiont", "free symbiont"], \
				 ["patient", "swab", "abscess", "lesion", "outbreak", "hospital", "clinical"], \
				 ["ATCC", "CIP", "NCTC", "collection"]]
				 
clinicaltags = [("blood",), ("cerebrospinal fluid",), ("mastitis",), ("groin",), ("tracheal", "throat"), ("sputum",), \
                ("diarrhea", "gastroenteritis"), ("respiratory", "lung"), ("urine", "urinary"), ("pleural",), ("baby", "neonatal")]
                
labtags = [["derived from", "variant of", "inactivated", "mutant", "mutated", "transposon", "laboratory"]]

ddatefields = [['day', 'month', 'year'], \
               ['year', 'month', 'day'], \
               ['month', 'year'], \
               ['year', 'month'], \
               ['year']]
ddatepatfields = dict(zip(ldatepats, ddatefields))

# pattern to separate assembly_id and assembly_name
reassgenbank = re.compile('(GC[AF]_[0-9]{9}\.[0-9])_(.+)')
reass = re.compile('(.+?\.[0-9])_(.+)')

## functions

def match_std_str(s, llsstdstr, default="other", rule='sin'):
	"""Given a list of equivalent trait names, return the standard name = the last of the list"""
	for lsstdstr in llsstdstr:
		for sstdstr in lsstdstr:
			if rule=='sin':
				if s in sstdstr: return lsstdstr[-1]
			elif rule=='ins':
				if sstdstr in s: return lsstdstr[-1]
			else:
				raise ValueError
	else:
		return default

def parse_assembly_name(assembname, reass=reass, group=0):
	seass = reass.match(assembname)
	geass = seass.groups()
	if group=='all': return geass
	else: return geass[group]
	
## script

opts, args = getopt.getopt(sys.argv[1:], '', ['assembly_folder_list=', 'add_raw_metadata=', 'add_curated_metadata=', 'add_dbxref=', 'add_assembly_info_dir=', 'default_species_name=', 'output='])
dopt = dict(opts)
nfldirassemb = dopt['--assembly_folder_list']
nfdhandmetaraw = dopt.get('--add_raw_metadata')
nfdhandmetacur = dopt.get('--add_curated_metadata')
nfdhanddbxref = dopt.get('--add_dbxref')
dirassemblyinfo = dopt.get('--add_assembly_info_dir')
defspename = dopt.get('--default_species_name')
output = dopt.get('--output')

## metadata extraction from compressed GenBank flat files
with open(nfldirassemb, 'r') as fldirassemb:
	ldirassemb = [line.rstrip('\n') for line in fldirassemb]


dmetadata = {}
lassembname = [os.path.basename(dirassemb) for dirassemb in ldirassemb]
lassemb = [parse_assembly_name(assembname, reass=reass) for assembname in lassembname]
lqualif = []

headerdbstart = "DBLINK      "
headerftstart = "FEATURES    "
headerpmstart = "   PUBMED   "
ftsourcestart = "     source "
smlemptystart = "            "
lrgemptystart = "                     "

for i, assembname in enumerate(lassembname):
	assemb = parse_assembly_name(assembname, reass=reass)
	qualif = None
	val = None
	with gzip.open("%s/%s_genomic.gbff.gz"%(ldirassemb[i], assembname), 'rb') as gbff:
		gbheader = True
		dbxrefblock = False
		sourceblock = 0
		lpmid = []
		ldbxref = []
		openpar = 0
		for line in gbff:
			if gbheader:   
				if line.startswith(headerftstart): 
					gbheader = False
					if lpmid: dmetadata.setdefault("pubmed_id", {})[assemb] = ','.join(lpmid)
					if ldbxref: dmetadata.setdefault("dbxref", {})[assemb] = ';'.join(ldbxref)
				if line.startswith(headerdbstart):
					dbxrefblock = True
				if line.startswith(headerpmstart):
					lpmid.append(line.strip('\n').split(headerpmstart)[-1])
				if dbxrefblock:
					for linestart in (headerdbstart, smlemptystart):
						if line.startswith(linestart):
							ldbxref.append(line.rstrip('\n').split(linestart)[1])
							break # for linestart loop
					else:
						dbxrefblock = False
			if sourceblock==1:
				# only captures the first source block ; second and more can be other organisms 
				# located within the bigger organism (e.g. inserted prophages) and should be ignored
				if line.startswith(lrgemptystart):
					li = line.strip()
					if li.startswith('/'):
						#~ print li
						qualval = li.strip('/\n').split('=', 1)
						if len(qualval)>1:
							# ignore qualifiers without value (e.g. '/focus' in accessions with multiple source blocks)
							qualif, val = qualval
							if not qualif in lqualif: lqualif.append(qualif) # Calife a la place du Calife!
							dmetadata.setdefault(qualif, {})[assemb] = val
					else:
						if qualif is None: continue # for line loop
						dmetadata[qualif][assemb] += ' '+li.strip('\n')
				else:
					break # for line loop
			if line.startswith(ftsourcestart):
				sourceblock += 1
	#~ print assemb, dmetadata['organism'][assemb]

# read additional data extracted by hand, if provided:
if nfdhandmetaraw:
	with open(nfdhandmetaraw, 'r') as fdhandmetaraw:
		for line in fdhandmetaraw:
			assemb, qualif, val = line.rstrip('\n').split('\t')
			dmetadata.setdefault(qualif, {})[assemb] = val
			
## read potential sequencing/assembly statistics and metadata
if dirassemblyinfo:
	lnfassemblyfiles = os.listdir(dirassemblyinfo)
	#~ for infotype in ['contig-N50', 'contig-count', 'assembly_level.tab', 'sequencing_technology.tab']:
	for infotype in ['assembly_level.tab', 'sequencing_technology.tab']:
		lnfinfo = [nf for nf in lnfassemblyfiles if nf.endswith(infotype)]
		infotag = infotype.replace('-', '_').replace('.tab', '')
		for nfinfo in lnfinfo:
			with open("%s/%s"%(dirassemblyinfo, nfinfo), 'r') as finfo:
				for line in finfo:
					assembname, val = line.rstrip('\n').split('\t')
					assemb = parse_assembly_name(assembname, reass=reass)
					dmetadata.setdefault(infotag, {})[assemb] = val.strip('" \r')

if os.path.exists(output):
        if os.path.isfile(output):
                raise ValueError, 'Output %s is an existing file'%output
else:
        os.mkdir(output)

nfout = os.path.join(output, 'metadata.tab')
with open(nfout, 'w') as fout:
	fout.write('\t'.join(['assembly']+lqualif)+'\n')
	for assemb in lassemb:
		fout.write('\t'.join([assemb]+[dmetadata[qualif].get(assemb, '') for qualif in lqualif])+'\n')

# separate many-to-many record of assemblies:dbxref, including pubmed_ids
nfoutdbxref = os.path.join(output, 'dbxrefs.tab')
with open(nfoutdbxref, 'w') as foutdbxref:
	for assemb in lassemb:
		if assemb in dmetadata["dbxref"]:
			dbxrefs = [dbxref.split(': ') for dbxref in dmetadata["dbxref"][assemb].split(';')]
			for db, xrefs in dbxrefs:
				if db!='Assembly': foutdbxref.write('\n'.join(['\t'.join([assemb, db, xref.strip()]) for xref in xrefs.split(',')])+'\n')
		if assemb in dmetadata["pubmed_id"]:
			pmids = dmetadata["pubmed_id"][assemb]
			foutdbxref.write('\n'.join(['\t'.join([assemb, "PubMed", pmid]) for pmid in pmids.split(',')])+'\n')
	## read file with manually added entries and copy-paste records
	if nfdhanddbxref:
		with open(nfdhanddbxref, 'r') as fhanddbxref:
			foutdbxref.write(fhanddbxref.read())

## curate the information
dcurated = {}
lcurqualif = ['organism', 'species', 'subspecies', 'serovar', 'strain', \
              'taxid', 'primary_pubmed_id', 'country', 'isolation_source', 'host', \
              'clinical_source', 'collection_year', 'collection_month', 'collection_day', 'sequencing_technology', \
              'sequencing_coverage', 'note']
              #~ 'sequencing_coverage', 'contig_N50', 'contig_count', 'note']
for cq in lcurqualif: dcurated[cq] = {}
na = ''
#~ na = 'NA'

for assemb in lassemb:
	coun = dmetadata.get('country',{}).get(assemb, na).split(':')[0]
	host = dmetadata.get('host',{}).get(assemb, na)
	isol = dmetadata.get('isolation_source',{}).get(assemb, na).lower().strip(' "')
	organism = dmetadata.get('organism',{}).get(assemb,na).strip(' "')
	strain = dmetadata.get('strain',{}).get(assemb, na)
	coldate = dmetadata.get('collection_date',{}).get(assemb, na)
	ppmid = dmetadata.get('pubmed_id',{}).get(assemb, na).split(',')[0]
	isolate = dmetadata.get('isolate',{}).get(assemb, na)
	sero = dmetadata.get('serovar',{}).get(assemb, na)
	note = dmetadata.get('note',{}).get(assemb, na)
	taxid = dict(dbxref.split(':') for dbxref in dmetadata.get('db_xref',{}).get(assemb,na).strip(' "').split(';'))['taxon']
	subspe = na
	ecol = na
	## strain
	# from organism name
	subspestr = organism.split()
	#~ print 'organism:', organism, 'subspestr:', subspestr,
	if 'symbiont' in organism:
		# specific case of organisms defined as symbionts of another
		ecol = 'symbiont'
		if subspestr[-1].endswith('symbiont'):
			# case 'Genus species symbiont' means it is the symbiont of Genus species, with no knowledge of its own species name
			host = ' '.join(subspestr[:-1])
			species = na
		elif 'symbiont of' in organism.lower():
			for symprefix in ['primary endo', 'secondary endo', 'endo', '']:
				symof = symprefix+'symbiont of'
				if symof in organism:
					species, host = organism.split(symof) ; host = host.strip(' ') ; species = species.strip(' ')
					parenmatch = re.match('^(.+) \(.+?\)$', host)
					if parenmatch:
						# expression in parenthesis in the end is more likely to refer to the symbion than the host,
						# e.g. in 'primary endosymbiont of Genus species (bacterium)' ; thus discard the parenthesis group
						host = host.group(0)
					break # for symprefix loop
	else:
		# generic prokaryote binomial nomenclature, with possible subspecies, serotype, strain extension
		species = ' '.join(subspestr[:2])
		subspestr = subspestr[2:]
		if len(subspestr)>1:
			if len(subspestr)>2 and subspestr[0]=='subsp.':
				subspe = subspestr[1]
				subspestr = subspestr[2:]
			if len(subspestr)>2 and subspestr[0]=='serovar' or subspestr[0]=="sv.":
				sero = subspestr[1]
				subspestr = subspestr[2:]
			if len(subspestr)>2 and subspestr[0]=='str.':
				strain = ' '.join(subspestr[1:])
			if strain==na:
				if isolate:
					strain = isolate
				else:
					strain = subspestr[-1]
	if species==na:
		species = defspename
	ssp = species.split()
	if (len(ssp)==1 or (len(ssp)==2 and ssp[0]=='Candidatus')) and ssp[-1][0].isupper():
		# name is just a 'Genus' name, add 'sp.'
		species = species+' sp.'
	#~ print 'ecol:', ecol, 'species:', species, 'host:', host
	dcurated['organism'][assemb] = organism
	dcurated['species'][assemb] = species
	dcurated['subspecies'][assemb] = subspe
	dcurated['serovar'][assemb] = sero
	dcurated['strain'][assemb] = strain
	dcurated['taxid'][assemb] = taxid
	## country of origin
	dcurated['country'][assemb] = coun
	## PubMed link
	dcurated['primary_pubmed_id'][assemb] = ppmid
	## host
	if host!=na:
		defhost = host if ecol=='symbiont' else 'other'
		host = match_std_str(host.lower(), \
								[["homo sapiens", "human"], \
								 ["bos taurus", "bovine"], \
								 ["mus musculus", "murine"], \
								 ["sus domesticus", "porcine"], \
								 ["equus", "horse"], \
								 ["gallus", "chicken"], \
								 ["canis", "dog"], \
								 ["musa", "banana"]], \
								 default=defhost)
	dcurated['host'][assemb] = host
	### isolation source
	## check occurence of clinical terms
	clinsource = set([])
	for tagtup in clinicaltags:
		for tag in tagtup:
			if tag in isol: clinsource.add(tagtup[0])
	clinical_source = ", ".join(list(clinsource))
	## filter out laboratory strains (i.e. experimental constructs or lab-adapted strains)
	for field in [ecol, isol, note]:
		if field and field!=na:
			ecol = match_std_str(field, labtags, default=ecol, rule='ins')
			if ecol=="laboratory":
				# strains to exclude!
				break # the for field loop
	else:
		# special case of E. coli lab/industrial strains
		if species=="Escherichia coli":
			for tag in ["K-12", "BL21", "DE3"]:
				if ((tag in strain) or (tag in isolate) or (tag in organism)):
					# strains to exclude! 
					ecol = "laboratory"
					break # the for tag loop
	# check in record notes potentially describing mutagenesis experiments 
	## match generic ecological source
	for field in [isol, note]:
		if ecol==na and field and field!=na:
			ecol = match_std_str(field, qualifgroups, default=na, rule='ins')
	dcurated['isolation_source'][assemb] = ecol
	dcurated['clinical_source'][assemb] = clinical_source
	# colection date
	dcurated['collection_year'][assemb] = na
	dcurated['collection_month'][assemb] = na
	dcurated['collection_day'][assemb] = na
	for datepat in ldatepats:
		res = datepat.search(coldate)
		if res:
			for i, g in enumerate(res.groups()):
				dcurated['collection_%s'%(ddatepatfields[datepat][i])][assemb] = dmonalias.get(g, g)
			break # for datepat	loop
	
	## sequencing info + note
	#~ for infotype in ['contig_N50', 'contig_count', 'sequencing_technology', 'note']:
	for infotype in ['sequencing_technology', 'note']:
		dcurated[infotype][assemb] = dmetadata[infotype].get(assemb, na)
	
	print '; '.join([assemb, species, strain, ppmid, ecol, note])
			
## read file with manually added entries and copy-paste records
if nfdhandmetacur:
	with open(nfdhandmetacur, 'r') as fmanualmetacur:
		for line in fmanualmetacur:
			assemb, qualif, val = line.rstrip('\n').split('\t')
			sval = val.strip(' "')
			if qualif=='note': 
				note = dcurated[qualif].setdefault(assemb, '')
				if note: dcurated[qualif][assemb] = (note+'; '+sval).replace("'", '').replace('"', '')
				else:dcurated[qualif][assemb] = sval.replace("'", '').replace('"', '')
			else: dcurated[qualif][assemb] = sval.replace("'", '').replace('"', '')

nfoutcur = os.path.join(output, 'metadata_curated.tab')
with open(nfoutcur, 'w') as foutcur:
	foutcur.write('\t'.join(['assembly_id', 'assembly_name']+lcurqualif)+'\n')
	for assembname in lassembname:
		geass = parse_assembly_name(assembname, reass=reass, group='all')
		assemb = geass[0]
		foutcur.write('\t'.join(list(geass)+[dcurated[qualif].get(assemb, '').strip('" ') for qualif in lcurqualif])+'\n')

