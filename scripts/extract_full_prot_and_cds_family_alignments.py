#!/usr/bin/python

import sys, getopt
import os
import subprocess
import multiprocessing
import gzip
import re
import time
import traceback
#~ from numpy import ndarray, zeros

#~ dryrun = True
dryrun = False

## functions
def nr2fullprotali(cdsfam, protid, dprotseq, dprotinfo, dprotfiletasks, dcdsfiletasks, dcdsfamsize, dcdstaxa):
	if not dryrun:
		# create full protein family alignment, i.e. with possibly redundant protein sequence (but unique id)
		fasta = dprotseq[protid]
		# extract protein product from fasta header
		patprod = re.compile(">%s (.+) \[.+\]"%(protid))
		sproduct = patprod.search(fasta[0])
		if not sproduct:
			patprod = re.compile(">%s (.+)$"%(protid))
			sproduct = patprod.search(fasta[0])
		if sproduct: product = " [product=%s]"%(sproduct.groups()[0].replace('MULTISPECIES: ', ''))
		else: product = ''
	# increment CDS count in family
	dcdsfamsize[cdsfam] = dcdsfamsize.get(cdsfam, 0) + len(dprotinfo[protid])
	ltasks = dprotinfo[protid]
	ltasks.sort()
	for task in ltasks:
		cdsid = task[1]
		# append extraction task
		# {src_cds_fasta_path:(unique_cds_id, cds_fam_id)}
		dcdsfiletasks.setdefault(task[0], []).append((cdsid, cdsfam))
		if not dryrun:
			# in fasta header, change nr prot id for unique cds_id, add taxon name and protein product
			fastahead = ">%s [%s]%s\n"%(cdsid, dcdstaxa[cdsid], product)
			dprotfiletasks[task] = [fastahead]+fasta[1:]
	#~ return ltasks

def loadsequences(nfnrprotaln):
	with open(nfnrprotaln, 'r') as fnrprotaln:
		protid = None
		dprotseq = {}
		for line in fnrprotaln:
			if line.startswith('>'):
				protid = line.strip('>\n').split()[0]
			dprotseq.setdefault(protid, []).append(line)
	return dprotseq
	
def sorttasksbysource(dprotfiletasks):
	# orderred tasks by source file, then by cds_id entry
	lprottasks = dprotfiletasks.keys()
	lprottasks.sort()
	return lprottasks
	
def joinlistasline(l):
	return '\t'.join([str(e) for e in l])+'\n'
	
def castPal2Nal(cdsfam, dirfullprotout, dirfullcdsseqout, fpal2nallog):
	"""sequential call to pal2nal.pl alignment reverse-translation script
	
	given gene family name, folders of protein alignements and aligned gene sequence, respectively,
	and a file handle to redirect pal2nal.pl verbose log; no return value"""
	nfprotali = "%s/%s.aln"%(dirfullprotout, cdsfam)
	nfcdsseq = "%s/%s.fasta"%(dirfullcdsseqout, cdsfam)
	p2ncmd = ["pal2nal.pl", "-output", "fasta", "-codontable", "11", nfprotali, nfcdsseq]
	#~ print ' '.join(p2ncmd)
	nfoutalnc = "%s/%s.aln"%(dirfullcdsaliout, cdsfam)
	with open(nfoutalnc, 'w') as foutalnc:
		subprocess.check_call(p2ncmd, stdout=foutalnc, stderr=fpal2nallog)
		
	
def castMultiPal2Nal(argtup):
	"""multi-threaded call to pal2nal.pl alignment reverse-translation script
	
	given a tuple containg the gene family name, folders of protein alignements and aligned gene sequence, respectively;
	returns pal2nal.pl verbose log (text string)"""
        try:
                # with multiprocessing
	        cdsfam, dirfullprotout, dirfullcdsseqout, fpal2nallog, queue = argtup
	        nfprotali = "%s/%s.aln"%(dirfullprotout, cdsfam)
	        nfcdsseq = "%s/%s.fasta"%(dirfullcdsseqout, cdsfam)
	        p2ncmd = ["pal2nal.pl", "-output", "fasta", "-codontable", "11", nfprotali, nfcdsseq]
	        #~ print ' '.join(p2ncmd)
	        nfoutalnc = "%s/%s.aln"%(dirfullcdsaliout, cdsfam)
	        foutalnc = open(nfoutalnc, 'w')
	        p2npipe = subprocess.Popen(p2ncmd, stdout=foutalnc, stderr=subprocess.PIPE)
	        foutalnc.close()
	        queue.put(cdsfam)
	        return p2npipe.stderr.read()
        except Exception, e:
                print "caught exception:"
                traceback.print_exc()
                sys.stdout.flush()
                raise e
	#~ fpal2nallog.write(p2npipe.stderr.read())
	
def castMultiBMGE(argtup):
	"""multi-threaded call to BMGE alignment trimming/filtering program
	
	Requires Java >1.6.
	"""
	pass

def main(dirnrprotaln, nfsingletonfasta, nfprotinfotab, nfreplinfotab, dirassemb, dirout, fam_prefix, dirlogs, nfidentseq=None):

	# define output folders
	suffdirout = os.path.basename(dirout)

	dirfullprotout = "%s/full_protfam_alignments"%(dirout)
	dirfullcdsseqout = "%s/full_cdsfam_fasta"%(dirout)
	dirfullcdsaliout = "%s/full_cdsfam_alignments"%(dirout)

	## non-redundant (nr) protein db
	prefixprotfam = fam_prefix + 'P'
	## (redundant) CDS db
	# diverse families
	prefixcdsfam = fam_prefix + 'C'
	#~ # homogeneous families (members of a family code a same unique protein) derived from ENTCGP000000; include singletons, i.e. ORFans pooled into ENTCGS000000
	#~ prefixsinglefam = 'ENTCGS'
	padlen = 6

	# extracted CDS sequences will be buffered in memory to avoid file open/close operations at every sequence
	# not to put too high as multiplied by the number of concurrently parsed gene families, i.e. potential the number of gene families
	nbseq2flushfam = 50

	# optional parsing of table of identical proteins
	didentseq = {}
	if nfidentseq:
		print "parse redundant protein names from '%s'"%nfidentseq
		curfam = None
		refprot = None
		with open(nfidentseq, 'r') as fidentseq:
			for line in fidentseq:
				fam, prot = line.rstrip('\n').split('\t')
				if fam == curfam:
					didentseq[prot] = refprot
				else:
					curfam = fam
					refprot = prot
	
	print "# parse replicon/genome assembly data"
	drepliasse = {}
	dreplitaxa = {}
	dassemtaxa = {}
	with open(nfreplinfotab, 'r') as freplinfotab:
		header = freplinfotab.readline()
		for line in freplinfotab:
			lsp = line.rstrip('\n').split('\t')
			# {refseq_acc:(assembly_id, assembly_name)}
			drepliasse[lsp[2]] = (lsp[0], lsp[1])
			# {refseq_acc:taxon_name}
			dreplitaxa[lsp[2]] = lsp[-1]
			# {assembly_id_name:(assembly_id, assembly_name, taxon_name)}
			dassemtaxa['_'.join(lsp[:2])] = tuple(lsp[:2]+lsp[-1:])

	dnfcdsfasta = {}
	dassemb = {}
	lnfcdsfasta = []

	print "# map genome assembly / CDS sequence dump files"
	for d in os.listdir(dirassemb):
		nfcdsfasta = os.path.join(dirassemb, d, d + "_cds_from_genomic.fna.gz")
		dnfcdsfasta[d] = nfcdsfasta
		dassemb[nfcdsfasta] = d
		lnfcdsfasta.append(nfcdsfasta)
	lnfcdsfasta.sort()

	print "# parse protein / CDS data"
	dprotinfo = {}
	#~ dcdsinfo = {}
	dcdstaxa = {}
	with open(nfprotinfotab, 'r') as fprotinfotab:
		header = fprotinfotab.readline()
		for line in fprotinfotab:
			lsp = line.rstrip('\n').split('\t')
			protid = lsp[0]
			repliacc = lsp[1]
			cdsid = lsp[-1]
			# {unique_cds_id:taxon_name <- refseq_acc}
			dcdstaxa[cdsid] = dreplitaxa[repliacc]
			# cds_fasta_path <- assembly_id_name <- refseq_acc
			nfcdsfasta = dnfcdsfasta['_'.join(drepliasse[repliacc])]
			#~ # {unique_cds_id:(cds_fasta_path, nr_prot_id)}
			#~ dcdsinfo[cdsid] = (nfcdsfasta, protid)
			# {nr_prot_id:[(cds_fasta_path, unique_cds_id), ...]}
			dprotinfo.setdefault(protid, []).append((nfcdsfasta, cdsid))
			
	if not os.path.exists(dirout):
		os.mkdir(dirout)
	if not os.path.exists(dirlogs):
		os.mkdir(dirlogs)
	if not os.path.exists(dirfullprotout):
		os.mkdir(dirfullprotout)
	if not os.path.exists(dirfullcdsseqout):
		os.mkdir(dirfullcdsseqout)
	if not os.path.exists(dirfullcdsaliout):
		os.mkdir(dirfullcdsaliout)

	# dict of lists of tasks of sequence extraction, orderred by source file
	dprotfiletasks = {}
	dcdsfiletasks = {}
	# parse nr protein alignments to generate protein alignment and task list for CDS sequence extration, to ensure order of sequences in both will be matching (required for pal2nal)
	lnfnrprotaln = ["%s/%s"%(dirnrprotaln, nf) for nf in os.listdir(dirnrprotaln)]
	# get last number in family id series for increment of new homogeneous CDSs families derived from the nr protein singleton one
	allprotfams = [os.path.basename(nf).rsplit('.', 1)[0] for nf in lnfnrprotaln]
	allprotfams.sort()
	lastfam = int(allprotfams[-1].split(prefixprotfam)[-1])
	homonfam = lastfam
	print "  last registered family id: %d"%lastfam

	dcdsfamsize = {}

	print "# parse singleton nr protein sequences"
	# first deal with the singleton nr proteins
	orfanfam = prefixcdsfam+(os.path.basename(nfsingletonfasta).rsplit('.', 1)[0].split(prefixprotfam)[-1])
	allcdsfam = []
	orfans = []
	dprotseq = loadsequences(nfsingletonfasta)
	dorfanprotfiletasks = {}
	for protid in dprotseq:
		# singleton nr protein
		if len(dprotinfo[protid])>1:
			# non-ORFan, homogeneous family
			homonfam += 1
			cdsfam = prefixcdsfam+(str(homonfam).zfill(padlen))
			allcdsfam.append(cdsfam)
			dprotfiletasks = {}
			# create full protein family alignment
			nr2fullprotali(cdsfam, protid, dprotseq, dprotinfo, dprotfiletasks, dcdsfiletasks, dcdsfamsize, dcdstaxa)
			if not dryrun:
				with open("%s/%s.aln"%(dirfullprotout, cdsfam), 'w') as foutfullprot:
					for prottask in sorttasksbysource(dprotfiletasks):
						foutfullprot.writelines(dprotfiletasks[prottask])
		else:
			# ORFan CDS
			# just add sequence to set of ORFans
			nr2fullprotali(orfanfam, protid, dprotseq, dprotinfo, dorfanprotfiletasks, dcdsfiletasks, dcdsfamsize, dcdstaxa)
			orfans.append(dprotinfo[protid][0][1]) # [unique_cds_id, ...]

	if not dryrun:
		with open("%s/%s.fasta"%(dirfullprotout, orfanfam), 'w') as forfanfullprot:
			for prottask in sorttasksbysource(dorfanprotfiletasks):
				forfanfullprot.writelines(dorfanprotfiletasks[prottask])	
		
	print "  %d final non-ORFan CDS families, including %d homogeneous (singleton protein derived) families"%(homonfam, homonfam-lastfam)
	print "  %d ORFan / %d total CDSs"%(dcdsfamsize[orfanfam], len(dcdstaxa))

	print "# parsing nr protein alignments"
	nnrprot = 0
	# then parse the regular diverse protein families
	for nfnrprotaln in lnfnrprotaln:
		protfam = os.path.basename(nfnrprotaln).rsplit('.', 1)[0]
		cdsfam = prefixcdsfam+(protfam.split(prefixprotfam)[-1])
		allcdsfam.append(cdsfam)
		dprotseq = loadsequences(nfnrprotaln)
		dprotfiletasks = {}
		# create full protein family alignment
		for protid in dprotseq:
			nr2fullprotali(cdsfam, protid, dprotseq, dprotinfo, dprotfiletasks, dcdsfiletasks, dcdsfamsize, dcdstaxa)
		if not dryrun:
			with open("%s/%s.aln"%(dirfullprotout, cdsfam), 'w') as foutfullprot:
				for prottask in sorttasksbysource(dprotfiletasks):
					foutfullprot.writelines(dprotfiletasks[prottask])				
		nnrprot += 1
		sys.stdout.write("\r%d / %d"%(nnrprot, len(lnfnrprotaln)))
		
	sys.stdout.write("\n")			

	lassemb = []
	print "# write gene family lists and genome gene content matrices"
	# write out genome list as used to order CDS and in gene content matrices
	with open("%s/genomes_ordered.tab"%(dirout), 'w') as foutgen:
		for nfcdsfasta in lnfcdsfasta:
			assembtax = dassemtaxa[dassemb[nfcdsfasta]]
			foutgen.write(joinlistasline(assembtax))
			lassemb.append(assembtax[0])

	allcdsfam.sort()
	# write down family contents
	# matrix genome x non-ORFan families
	dfamassembcount = dict([(cdsfam, [0]*len(lnfcdsfasta)) for cdsfam in allcdsfam])
	#~ afamassembcount = zeros((len(lnfcdsfasta), len(allcdsfam)))
	# matrix genome x ORFans (unique CDS ids)
	dorfanassembcount = dict([(orfan, [0]*len(lnfcdsfasta)) for orfan in orfans])

	foutfamcontent = open("%s/full_families_info-noORFans.tab"%(dirout), 'w')
	foutorfancontent = open("%s/%s_info-ORFans.tab"%(dirout, orfanfam), 'w')
	for i, nfcdsfasta in enumerate(lnfcdsfasta):
		assemblyidname = dassemb[nfcdsfasta]
		assemblyid = '_'.join(assemblyidname.split('_', 2)[:2])
		for task in dcdsfiletasks[nfcdsfasta]:
			cdsid, cdsfam = task
			if cdsfam==orfanfam:
				foutorfancontent.write(joinlistasline([assemblyid, cdsid, cdsfam]))
				dorfanassembcount[cdsid][i] += 1
			else:
				foutfamcontent.write(joinlistasline([assemblyid, cdsid, cdsfam]))
				dfamassembcount[cdsfam][i] += 1
	foutfamcontent.close()
	foutorfancontent.close()
	print "matrix of family counts( genome x non-ORFan families):"
	nfoutfammat = "%s/full_families_genome_counts-noORFans.mat"%(dirout)
	print "  '%s'"%nfoutfammat
	with open(nfoutfammat, 'w') as foutfammat:
		foutfammat.write(joinlistasline(['']+lassemb))
		for cdsfam in allcdsfam:
			foutfammat.write(joinlistasline([cdsfam]+dfamassembcount[cdsfam]))

	print "matrix of family counts( genome x ORFan families):"
	nfoutorfanmat = "%s/%s_genome_counts-ORFans.mat"%(dirout, orfanfam)
	print "  '%s'"%nfoutorfanmat
	with open(nfoutorfanmat, 'w') as foutorfanmat:
		foutorfanmat.write(joinlistasline(['']+lassemb))
		for orfan in orfans:
			foutorfanmat.write(joinlistasline([orfan]+dorfanassembcount[orfan]))
		 
	# search for core families
	with open("%s/core-genome_families.tab"%(dirout), 'w') as foutcore:
		for cdsfam in allcdsfam:
			famcounts = dfamassembcount[cdsfam]
			for i in range(len(famcounts)):
				if famcounts[i]!=1:
					#~ print i, "(%d)"%famcounts[i],
					break
			else:
				# core family
				foutcore.write(cdsfam+'\n')

	# tasks are sorted by source file, must then order them by cds_id entry
	for nfcdsfasta in dcdsfiletasks:
		dcdsfiletasks[nfcdsfasta].sort(key=lambda x: x[0])
	# to avoid file open/close operations at every sequence, without having to maintain many open connections, buffer content in memory
	dfambuffer = {}
	# keep track of how many sequences already extracted for a family so can flush buffer when reaching threshold or done
	dnextractcds = {}
	ncdssrc = 0

	if not dryrun:
		print "# extract CDSs from genomic dump files"
		## proceed by source file to extract CDSs
		# the files are ordered by their path, as should be the CDS entries in each family
		for nfcdsfasta in lnfcdsfasta:
			with gzip.open(nfcdsfasta, 'rb') as fcdsfasta:
				cdsid = None
				dcdsseq = {}
				for line in fcdsfasta:
					if line.startswith('>'):
						cdsid = line.strip('>\n').split()[0]
					dcdsseq.setdefault(cdsid, []).append(line)
			ncdssrc += 1
			#~ print nfcdsfasta
			for task in dcdsfiletasks[nfcdsfasta]:
				cdsid, cdsfam = task
				#~ print '', cdsid, cdsfam
				dfambuffer.setdefault(cdsfam, []).append(''.join(dcdsseq[cdsid]))
				#~ dfambuffer.setdefault(cdsfam, {}).update({cdsid:''.join(dcdsseq[cdsid])})
				nprevextractcds = dnextractcds.get(cdsfam, 0)
				nextrbuffcds = nprevextractcds + len(dfambuffer[cdsfam])
				if len(dfambuffer[cdsfam])>=nbseq2flushfam or nextrbuffcds==dcdsfamsize[cdsfam]:
					# all CDS extracted, flush buffer
					if nprevextractcds>0: wa = 'a'
					else: wa = 'w'
					with open("%s/%s.fasta"%(dirfullcdsseqout, cdsfam), wa) as fdest:
						fdest.writelines(dfambuffer[cdsfam])
					# increment extracted CDS count in family
					dnextractcds[cdsfam] = nextrbuffcds
					del dfambuffer[cdsfam]
			sys.stdout.write("\r%d\tsource files parsed ; %d\tfamilies in buffer"%(ncdssrc, len(dfambuffer)))

		sys.stdout.write("\n")
		if dfambuffer:
			# all connections should be shut at the end of the loop
			raise ValueError, "not all destination CDS family fasta file connections are closed ; some some family must expect more sequences"
		
		## cast pal2nal on (full protein alignment, unaligned CDS fasta) file pairs
		# detect the number of processors
		print "# reverse translate alignments"
		dirpal2nallogs = "%s/pal2nal"%dirlogs
		if not os.path.exists(dirpal2nallogs):
			os.mkdir(dirpal2nallogs)
		fpal2nallog = open("%s/extract_full_prot+cds_family_alignments.pal2nal_log"%dirpal2nallogs, 'w')
		# sequentially
		#~ np2p = 0
		#~ for cdsfam in allcdsfam:
			#~ castPal2Nal(cdsfam, dirfullprotout, dirfullcdsseqout, fpal2nallog)
			#~ np2p += 1
			#~ sys.stdout.write("\r%d\t/%d"%(np2p, len(allcdsfam)))
		
		# with multithreading ; using map_assync and a queue, untested
		nbcores = multiprocessing.cpu_count()
		pool = multiprocessing.Pool(processes=nbcores)
		manag = multiprocessing.Manager()
		queue = manag.Queue()
		p2nlog = pool.map(castMultiPal2Nal, ((cdsfam, dirfullprotout, dirfullcdsseqout, fpal2nallog, queue) for cdsfam in allcdsfam))
		#~ result = pool.map_async(castMultiPal2Nal, ((cdsfam, dirfullprotout, dirfullcdsseqout, fpal2nallog, queue) for cdsfam in allcdsfam))
		#~ # monitor loop
		#~ while True:
			#~ if result.ready():
				#~ break
			#~ else:
				#~ sys.stdout.write("\r%d\t/%d"%(queue.qsize(), len(allcdsfam)))
				#~ time.sleep(0.5)

		#~ p2nlog = result.get()
		fpal2nallog.write('\n'.join(p2nlog)+'\n')

		fpal2nallog.close()
		sys.stdout.write("\n")

def usage():
	modulename = 'extract_full_prot_and_cds_family_alignments.py'
	try:
		nfscript = sys.argv[0] if modulename in sys.argv[0] else modulename
	except IndexError:
		nfscript = modulename
	s =  'Usage:\n'
	s += 'python %s '%nfscript
	s += 'Mandatory options:\n'
	s += '--nrprot_fam_alns\t\tpath to folder containing all the non-redundant protein family alignments (FASTA format)\n'
	s += '--singletons\t\tpath to file of __unaligned__ singleton protein sequences (FASTA format)\n'
	s += '--prot_info\t\tpath to file of all protein/CDS informations (table as generated by pantagruel/scripts/allgenome_gff2db.py)\n'
	s += '--repli_info\t\tpath to file of all assembly/replicon_accession informations (table as generated by pantagruel/scripts/allgenome_gff2db.py)\n'
	s += '--assemblies\t\tpath to folder containing all assemblies (GenBank FTP format and file structure)\n'
	s += '--dirout\t\tpath to output folder\n'
	s += '--famprefix\t\tstring used sa prefix for protein and CDS families (without the trailing \'C\' or \'P\' and family id numerals)\n'
	s += 'Other options:\n'
	s += '--logs\t\tpath to log folder (default to \'dirout/logs\')\n'
	s += '--identical_seqs\t\tpath to table of identical protein sequences'
	return s




## script
if __name__ == '__main__':

	opts, args = getopt.getopt(sys.argv[1:], 'h', ['nrprot_fam_alns=', 'assemblies=', 'singletons=', \
	                                               'prot_info=', 'repli_info=', 'identical_prots=', \
	                                               'famprefix=', \
	                                               'dirout=', 'logs=', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	dirnrprotaln = dopt['--nrprot_fam_alns'].rstrip('/')
	nfsingletonfasta = dopt['--singletons']
	nfprotinfotab = dopt['--prot_info']
	nfreplinfotab = dopt['--repli_info']
	dirassemb = dopt['--assemblies']
	dirout = dopt['--dirout']
	fam_prefix = dopt['--famprefix']
	dirlogs = dopt.get('--logs', "%s/logs"%(dirout))
	nfidentseq = dopt.get('--identical_prots')
	
	main(dirnrprotaln, nfsingletonfasta, nfprotinfotab, nfreplinfotab, dirassemb, dirout, fam_prefix, dirlogs, nfidentseq)
	
        
       
