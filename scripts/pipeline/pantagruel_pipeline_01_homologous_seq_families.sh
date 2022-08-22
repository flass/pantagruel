#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

checkptgversion
checkfoldersafe ${seqdb}

if [ ! -z "${ptgthreads}" ] ; then
  mmthreads="--threads ${ptgthreads}"
else
  mmthreads=""
fi

#############################
## 01. Homologous Sequence db
#############################
if [ -z "${ptgtmp}" ] ; then
  echo "Error: \${ptgtmp} is not defined; exit now"
  exit 1
fi
mmseqstmp=${ptgtmp}/mmseqs && rm -rf ${mmseqstmp} && mkdir -p ${mmseqstmp}/

if [ ! -z "${updatedbfrom}" ] ; then
  echo "this patagruel db will be generated by an updating a previously built one, located at '${updatedbfrom}'"
  prevdbseqdb=${updatedbfrom}/$(basename ${seqdb})
  if [ ! -d "${prevdbseqdb}" ] ; then
    echo "Error: could not find the protein clustering folder (output of Pantagruel task 01) from the reference Pantagruel db/update source where expected:"
	ls ${prevdbseqdb}
	exit 1
  else
    prevdballfaarad=${prevdbseqdb}/$(basename ${allfaarad})
	if [[ ! -e "${prevdballfaarad}.faa" || ! -e "${prevdballfaarad}.mmseqsdb" ]] ; then
      echo "Error: could not find the concatenated proteome fasta file and/or MMSeqs protein db from the reference Pantagruel db/update source where expected:"
	  ls ${prevdballfaarad}.faa ${prevdballfaarad}.mmseqsdb
	  exit 1
    fi
  fi
  mmsdbtag=".updated"
else
  mmsdbtag=""
fi

## extract all the protein sequences into single proteome fasta files
if [ -z "${allfaarad}" ] ; then
  echo "Error: \${allfaarad} is not defined; exit now"
  exit 1
fi
rm -f ${allfaarad}*
for ass in `ls ${assemblies}` ; do
 if [ -d ${assemblies}/${ass} ] ; then
  faa=$(ls ${assemblies}/${ass}/*protein.faa 2> /dev/null)
  if [[ ! -z "${faa}" && -s "${faa}" ]] ; then
   cat ${faa} >> ${allfaarad}.faa && echo ${faa} >> ${allfaarad}_list
  else
   faagz=$(ls ${assemblies}/${ass}/*protein.faa.gz 2> /dev/null)
   if [ -z "${faagz}" ] ; then
     echo "Error: could not find any proteome file matching '${assemblies}/${ass}/*protein.faa[.gz]' ; exit now"
     exit 1
   fi
   zcat ${faagz} >> ${allfaarad}.faa && echo ${faagz} >> ${allfaarad}_list
  fi
 fi
done
promptdate "-- $(wc -l ${allfaarad}_list | cut -d' ' -f1) proteomes in dataset"
promptdate "-- $(grep -c '>' ${allfaarad}.faa) proteins in dataset"

# dereplicate proteins in db based on their identifier as can occur among genomes from the same source RefSeq or Genbank 
python2.7 ${ptgscripts}/dereplicate_fasta.py ${allfaarad}.faa ${allfaarad}.nrprotids.faa
promptdate "-- $(grep -c '>' ${allfaarad}.nrprotids.faa) non-redundant protein ids in dataset"
diff ${allfaarad}.faa ${allfaarad}.nrprotids.faa > ${allfaarad}.fullVSnrprotids.faa.diff
if [ -z ${allfaarad}.fullVSnrprotids.faa.diff ] ; then
  echo "could not find same-name protein redundancy in full proteome file '${allfaarad}'; delete the '${allfaarad}.nrprotids.faa' to save disk space and instead make a symlink to full proteome '${allfaarad}'"
  rm ${allfaarad}.nrprotids.faa
  ln -s ${allfaarad}.faa ${allfaarad}.nrprotids.faa
fi
rm ${allfaarad}.fullVSnrprotids.faa.diff

mmseqslogs=${ptglogs}/mmseqs && mkdir -p ${mmseqslogs}/
## clustering of identical protein sequences
# notably those from the custom assemblies to those from the public database (and those redudant between RefSeq, Genbank or custom sets)

# run mmseqs clusthash with 100% seq id threshold
echo "${datepad}-- Perform first protein clustering step (100% prot identity clustering with clusthash algorithm)"
mmlog0=${mmseqslogs}/mmseqs-0-identicalprot-clusthash.log
mmseqs createdb ${allfaarad}.nrprotids.faa ${allfaarad}.mmseqsdb &> ${mmlog0}
if [ -z "${updatedbfrom}" ] ; then
	mmseqs clusthash ${mmthreads} --min-seq-id 1.0 ${allfaarad}.mmseqsdb ${allfaarad}.clusthashdb_minseqid100 &>> ${mmlog0}
	mmseqs clust ${mmthreads} ${allfaarad}.mmseqsdb ${allfaarad}.clusthashdb_minseqid100 ${allfaarad}.clusthashdb_minseqid100_clust &>> ${mmlog0}
else
	# update previous clustering
	mmseqs clusterupdate ${mmthreads} --min-seq-id 1.0 ${prevdballfaarad}.mmseqsdb ${allfaarad}.mmseqsdb ${prevdballfaarad}.clusthashdb_minseqid100_clust ${allfaarad}.mmseqsdb${mmsdbtag} ${allfaarad}.clusthashdb_minseqid100_clust ${mmseqstmp} &>> ${mmlog0}
fi
mmsummary0=$(tail -n 4 ${mmlog0} | head -n 3)
mmseqs createseqfiledb ${mmthreads} ${allfaarad}.mmseqsdb${mmsdbtag} ${allfaarad}.clusthashdb_minseqid100_clust ${allfaarad}.clusthashdb_minseqid100_clusters &>> ${mmlog0}
checkexec "First protein clustering step failed; please investigate error reports in '${mmlog0}'" "${datepad}-- First protein clustering step complete: ${mmsummary0}"

# get table of redundant protein names
if [ -z "${updatedbfrom}" ] ; then
    python2.7 ${ptgscripts}/split_mmseqs_clustdb_fasta.py ${allfaarad}.clusthashdb_minseqid100_clusters "NRPROT" ${allfaarad}.clusthashdb_minseqid100_families 6 0 0
else
    python2.7 ${ptgscripts}/split_mmseqs_clustdb_fasta.py ${allfaarad}.clusthashdb_minseqid100_clusters "NRPROT" ${allfaarad}.clusthashdb_minseqid100_families 6 0 0 ${prevdballfaarad}.identicals.tab
fi
grep -v NRPROT000000 ${allfaarad}.clusthashdb_minseqid100_families.tab > ${allfaarad}.identicals.tab
python2.7 ${ptgscripts}/genefam_table_as_list.py ${allfaarad}.identicals.tab ${allfaarad}.identicals.list 0
if [ -z "${updatedbfrom}" ] ; then
  python2.7 ${ptgscripts}/remove_identical_seqs.py ${allfaarad}.nrprotids.faa ${allfaarad}.identicals.tab ${allfaarad}.nr.faa
else
  python2.7 ${ptgscripts}/remove_identical_seqs.py ${allfaarad}.nrprotids.faa ${allfaarad}.identicals.tab ${allfaarad}.nr.faa ${prevdballfaarad}.identicals.tab
fi
## collect data from assemblies, including matching of (nr) protein to CDS sequence ids
python2.7 ${ptgscripts}/allgenome_gff2db.py --assemb_list ${genomeinfo}/assemblies_list --dirout ${genomeinfo}/assembly_info \
 --ncbi_taxonomy ${ncbitax} --identical_prots ${allfaarad}.identicals.list

## check consistency of non-redundant protein sets
mkdir -p ${ptgtmp}
protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'nr_protein_id' | cut -d':' -f1)
if [ -z "${protidfield}" ] ; then 
 protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'protein_id' | cut -d':' -f1)
fi
cut -f ${protidfield} ${genomeinfo}/assembly_info/allproteins_info.tab | grep -v "^$\|protein_id" | sort -u > ${genomeinfo}/assembly_info/allproteins_in_gff
grep '>' ${allfaarad}.nr.faa | cut -d' ' -f1 | cut -d'>' -f2 | sort -u > ${allfaarad}.nr_protlist
# compare original dataset of nr protein (as described in the input GFF files) to the aligned nr proteome
diff ${genomeinfo}/assembly_info/allproteins_in_gff ${allfaarad}.nr_protlist > ${ptgtmp}/diff_prot_info_fasta
if [ -s ${ptgtmp}/diff_prot_info_fasta ] ; then 
  >&2 promptdate
  >&2 echo "ERROR: inconsistent propagation of the protein dataset:"
  >&2 echo "present in combined proteome fasta files / absent in info table generated from input GFF:"
  >&2 grep '>' ${ptgtmp}/diff_prot_info_fasta | cut -d' ' -f2
  >&2 echo "present in info table generated from input GFF / absent in combined proteome fasta files:"
  >&2 grep '<' ${ptgtmp}/diff_prot_info_fasta | cut -d' ' -f2
  exit 1
fi

## clustering of proteome db with  MMSeqs2 
# (https://github.com/soedinglab/MMseqs2,  Steinegger M and Soeding J. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv, doi: 10.1101/079681 (2016))
# compute the memory use of MMSeq2: M = (7 × N × L + 8 × a^k) bytes, N the number of sequences, L their average size, a the size of the alphabet
## clustering of nr proteome 
# run mmseqs cluster with default parameters
# used MMseqs2 Version: e5d64b2701789e7eef8fcec0812ccb910c8dfef3
# compute the memory use of MMSeq2: M = (7 × N × L + 8 × a^k) bytes, N the number of sequences, L their average size, a the size of the alphabet
# create MMseqs2 db
mmlog1=${mmseqslogs}/mmseqs-1-cluster.log
mmseqs createdb ${allfaarad}.nr.faa ${allfaarad}.nr.mmseqsdb${mmsdbtag} &> ${mmlog1}
# perform clustering
mkdir -p ${families}
echo "${datepad}-- Perform second protein clustering step (to find homologs with cluster algorithm)"
mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
# perform similarity search and clustering
if [ -z "${updatedbfrom}" ] ; then
	mmseqs cluster ${mmthreads} ${allfaarad}.nr.mmseqsdb${mmsdbtag} $mmseqsclout $mmseqstmp &>> ${mmlog1}
else
	# update previous clustering
	prevdbmmseqsclout=${prevdbseqdb}/$(basename ${families})/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
	mmseqs clusterupdate ${mmthreads} ${prevdballfaarad}.nr.mmseqsdb ${allfaarad}.nr.mmseqsdb $prevdbmmseqsclout ${allfaarad}.nr.mmseqsdb${mmsdbtag} $mmseqsclout $mmseqstmp &>> ${mmlog1}
fi
mmsummary=$(tail -n 4 ${mmlog1} | head -n 3)
# generate indexed fasta file listing all protein families
mmseqs createseqfiledb ${mmthreads} ${allfaarad}.nr.mmseqsdb${mmsdbtag} $mmseqsclout ${mmseqsclout}_clusters &>> ${mmlog1}
checkexec "Second protein clustering step failed; please investigate error reports in '${mmlog1}'" "${datepad}-- Second protein clustering step complete: ${mmsummary1}"
# generate separate fasta files with family identifiers distinc from representative sequence name
python2.7 ${ptgscripts}/split_mmseqs_clustdb_fasta.py ${mmseqsclout}_clusters "${famprefix}P" ${mmseqsclout}_clusters_fasta 6 1 0
checkexec "Failed to split mmseqs cluster '${mmseqsclout}_clusters'" "${datepad}-- Successfully split mmseqs cluster '${mmseqsclout}_clusters'"
promptdate "-- $(wc -l ${mmseqsclout}_clusters_fasta.tab | cut -d' ' -f1) non-redundant proteins"
promptdate "-- classified into $(ls ${mmseqsclout}_clusters_fasta/ | wc -l) clusters"
echo "${datepad}-- including artificial cluster ${famprefix}P000000 gathering $(grep -c '>' ${mmseqsclout}_clusters_fasta/${famprefix}P000000.fasta) ORFan nr proteins"
echo "${datepad}-- (NB: some are not true ORFans as can be be present as identical sequences in several genomes)"

rm -rf ${mmseqstmp}
