#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

export ptgroot=$1
envsourcescript=${ptgroot}/environ_pantagruel.sh
source $envsourcescript

#############################
## 01. Homologous Sequence db
#############################

## extract all the protein sequences into single proteome fasta files
mkdir -p ${seqdb}
faafiletag='protein'
#~ faafiletag='translated_from_cds'

### TO DO ###
# using translated_from_cds could probably save some steps, 
# even though would mean handling slightly bigger files at the first clustering step (clusthash_minseqid100)
# after that, protein sequence files are already non-redundant step anyway
# also would ensure consistency of protein and CDS sequence at the pal2nal step

#~ export allfaarad=${seqdb}/all_proteomes
#~ rm -f ${allfaarad}*
#~ for ass in `ls ${assemblies}` ; do
 #~ faa=$(ls ${assemblies}/${ass}/ | grep '${faafiletag}.faa')
 #~ zcat $faa >> ${allfaarad}.faa && echo $faa >> ${allfaarad}_list
#~ done
#~ wc -l ${allfaarad}_list

ls ${indata}/assemblies/GC[AF]_* -d > ${genomeinfo}/refseq_assemblies_list
sed -e 's#\(.\+\)/\([^/]\+$\)#\1/\2/\2_protein\.faa\.gz#g' ${genomeinfo}/refseq_assemblies_list > ${faarefseq}_list
for faa in `cat ${faarefseq}_list` ; do zcat $faa >> ${faarefseq} ; echo $faa ; done 
grep -c '>' ${faarefseq}
# dereplicate proteins in db
nrfaarefseq=$(python $dbscripts/dereplicate_fasta.py $faarefseq)
echo "$(dateprompt)-- $(grep -c '>' $nrfaarefseq) non-redundant proteins in dataset"

export allfaarad=${seqdb}/all_proteomes
if [[ "$(ls -A "${customassemb}/contigs/" 2>/dev/null)" ]] ; then
  # combine Refseq and custom protomes
  cp -f $nrfaarefseq ${allfaarad}.faa
  for ass in `ls ${annot}/` ; do
    cat ${annot}/$ass/*faa >> ${allfaarad}.faa
  done
else
  mv $nrfaarefseq ${allfaarad}.faa
fi
grep -c '>' ${allfaarad}.faa

## clustering of identical protein sequences
# notably those from the custom assemblies to those from the public database (and those redudant between RefSeq and Genbank sets)
# run mmseqs clusthash with 100% seq id threshold
# used MMseqs2 Version: 6306925fa9ae6198116c26e605277132deff70d0
mmseqs createdb ${allfaarad}.faa  ${allfaarad}.mmseqsdb
mmseqs clusthash --min-seq-id 1.0 ${allfaarad}.mmseqsdb ${allfaarad}.clusthashdb_minseqid100
mmseqs clust ${allfaarad}.mmseqsdb ${allfaarad}.clusthashdb_minseqid100 ${allfaarad}.clusthashdb_minseqid100_clust
mmseqs createseqfiledb ${allfaarad}.mmseqsdb ${allfaarad}.clusthashdb_minseqid100_clust ${allfaarad}.clusthashdb_minseqid100_clusters

# get table of redundant protein names
python ${ptgscripts}/split_mmseqs_clustdb_fasta.py ${allfaarad}.clusthashdb_minseqid100_clusters "NRPROT" ${allfaarad}.clusthashdb_minseqid100_families 6 0 0
grep -v NRPROT000000 ${allfaarad}.clusthashdb_minseqid100_families.tab > ${allfaarad}.identicals.tab
python ${ptgscripts}/genefam_table_as_list.py ${allfaarad}.identicals.tab ${allfaarad}.identicals.list 0
python ${ptgscripts}/remove_identical_seqs.py ${allfaarad}.faa ${allfaarad}.identicals.list ${allfaarad}.nr.faa

## collect data from assemblies, including matching of (nr) protein to CDS sequence ids
python ${ptgscripts}/allgenome_gff2db.py --assemb_list ${genomeinfo}/assemblies_list --dirout ${genomeinfo}/assembly_info \
 --ncbi_taxonomy ${ncbitax} --identical_prots ${allfaarad}.identicals.list

## clustering of proteome db with  MMSeqs2 
# (https://github.com/soedinglab/MMseqs2,  Steinegger M and Soeding J. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv, doi: 10.1101/079681 (2016))
# compute the memory use of MMSeq2: M = (7 × N × L + 8 × a^k) bytes, N the number of sequences, L their average size, a the size of the alphabet
## clustering of nr proteome 
# run mmseqs cluster with default parameters
# used MMseqs2 Version: e5d64b2701789e7eef8fcec0812ccb910c8dfef3
# compute the memory use of MMSeq2: M = (7 × N × L + 8 × a^k) bytes, N the number of sequences, L their average size, a the size of the alphabet
mmseqslogs=${ptglogs}/mmseqs && mkdir -p $mmseqslogs
# create MMseqs2 db
mmseqs createdb ${allfaarad}.nr.faa ${allfaarad}.nr.mmseqsdb &> $mmseqslogs/mmseqs-createdb.log
# perform clustering
mmseqstmp=${ptgtmp}/mmseqs && rm -rf $mmseqstmp && mkdir -p $mmseqstmp
export families=${seqdb}/protein_families && mkdir -p $families
mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
# perform similarity search and clustering ; uses all CPU cores by default
mmseqs cluster ${allfaarad}.nr.mmseqsdb $mmseqsclout $mmseqstmp &> $mmseqslogs/$(basename $mmseqsclout).log
# generate indexed fasta file listing all protein families
mmseqs createseqfiledb ${allfaarad}.nr.mmseqsdb $mmseqsclout ${mmseqsclout}_clusters
# generate separate fasta files with family identifiers distinc from representative sequence name
python ${ptgscripts}/split_mmseqs_clustdb_fasta.py ${mmseqsclout}_clusters "${famprefix}P" ${mmseqsclout}_clusters_fasta 6 1 0
echo "$(dateprompt)-- $(wc -l ${mmseqsclout}_clusters_fasta.tab | cut -d' ' -f1) non-redundant proteins"
echo "$(dateprompt)-- classified into $(ls ${mmseqsclout}_clusters_fasta/ | wc -l) clusters"
echo "${datepad}-- including artificial cluster ${famprefix}P000000 gathering $(grep -c '>' ${mmseqsclout}_clusters_fasta/${famprefix}P000000.fasta) ORFan nr proteins"
echo "${datepad}-- (NB: some are not true ORFans as can be be present as identical sequences in several genomes)"
