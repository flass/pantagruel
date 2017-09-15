#!/bin/bash

#########################################################
## REC-A-PANG: 
##             a pipeline for 
##             phylogenetic reconciliation
##             of a bacterial pangenome
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 22 August 2017

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

#### Set environment variables
export myemail='me.myself@respectable-institu.ti.on'
export raproot='/path/to/base/folder/for/all/this/business'
export dbscripts='/path/to/pipeline/scripts'

export rapdbname='aBoldDatabaseName' # mostly name of the top folder
export famprefix='ABCDE'             # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.
# generic env var derived from the above
source ~/environ_rec-a-pang.sh

################################
## 00. Data download and grooming
################################

### Download data
## using FTP
# downloading relevant tables from NCBI Taxonomy database 
user='anonymous'
pswd=$myemail
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
source="/pub/taxonomy"
files="taxcat.tar.gz* taxcat_readme.txt taxdump.tar.gz* taxdump_readme.txt"
ncbitax="/enterobac/NCBI/Taxonomy"
mkdir -p $ncbitax
lftp -c "open $openparam; cd $source ; mget -O $ncbitax/ $files ; quit"
mv $ncbitax/README $ncbitax/readme_accession2taxid
# reduce taxonomic id2name reference file complexity
for tgz in `ls $ncbitax/*.tar.gz` ; do md5sum -c $tgz.md5 && tar -xzf $tgz ; done
grep "scientific name"  $ncbitax/names.dmp | sed -e 's/\t|\t/\t/g' | cut -f1,2 > $ncbitax/scientific_names.dmp

## using NCBI web interface:
# Search Assembly database using a query defined as convenient:	
# e.g.: 'txid543[Organism:exp] AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter]) AND ("complete genome"[filter])'
# save assemblies (Download Assemblies > Source Database: RefSeq; File type: All file type (including assembly-structure directory)) as
genome_assemblies.tar
# extract assembly data
assfolder=`tar -tf genome_assemblies.tar | grep ncbi-genomes | head -n 1`
tar -C $assemblies/ -xf genome_assemblies.tar ncbi-genomes-*/*
rm genome_assemblies.tar 
# store in centralised folder for NCBI assemblies and just record links
mv $indata/$assfolder/* $ncbiass/
mv $indata/$assfolder/ $indata/assemblies/
for ass in `cut -f3 $indata/assembly_result.txt | tail -n +2` ; do ln -s $ncbiass/${ass}_* $indata/assemblies/  ; done


### Groom data
## generate assembly statistics to verify genome finishing status
cd ${assemblies}/
mkdir -p ${indata}/assembly_stats/
filetag="_assembly_stats.txt"
patasslev="Assembly level"
patseqteq="Sequencing technology"
patcontig="all\tall\tall\tall\tcontig-"	
# assembly status
sedfilter0="s@\(.\+\)/.\+${filetag}:# ${patasslev}: \(.\+\)@\1\t\2@g"
grep "$patasslev" */*${filetag} | sed -e "$sedfilter0" | sed -e 's/\r//g' > ${indata}/assembly_stats/assembly_level.tab

mkdir ${indata}/complete_genomes/
rm -f ${allcomplete}*
ls ${assemblies}/GCF_* -d > ${allcomplete}_withlabstrains_list
# complete genomes do not have their sequencing / assembly details reported in *_assembly_stats.txt files ; rather lok for the sequencing metadata in the GenBank flat file
for assemb in `cat ${allcomplete}_withlabstrains_list` ; do 
ass=$(basename $assemb) ; tech=`zcat $assemb/${ass}_genomic.gbff.gz | grep -m 1 "Sequencing Technology" | awk -F' :: ' '{print $NF}' | sed -e 's/Il;lumina/Illumina/g'`
if [ -z "$tech" ] ; then tech='NA' ; fi
echo -e "${ass}\t${tech}" ; done > ${indata}/sequencing_technology.tab

## extract assembly/sample metadata from flat files
python $dbscripts/extract_metadata_from_gbff.py --assembly_folder_list=${allcomplete}_withlabstrains_list --add_raw_metadata=${complete}/manual_metadata_dictionary.tab \
--add_curated_metadata=${complete}/manual_curated_metadata_dictionary.tab --add_dbxref=${complete}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats
# remove lab strains
grep "laboratory" ${allcomplete}_withlabstrains_metadata_curated.tab | cut -f1 > $complete/labstrains_list
for strain in `cat $complete/labstrains_list` ; do
 grep $strain ${allcomplete}_withlabstrains_list >> ${allcomplete}_labstrains_list
done
diff ${allcomplete}_withlabstrains_list ${allcomplete}_labstrains_list | grep '<' | cut -d' ' -f2 > ${allcomplete}_list
# regenerate metadata tables without the lab strains
python $dbscripts/extract_metadata_from_gbff.py --assembly_folder_list=${allcomplete}_list --add_raw_metadata=${complete}/manual_metadata_dictionary.tab \
--add_curated_metadata=${complete}/manual_curated_metadata_dictionary.tab --add_dbxref=${complete}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats

#############################
## 01. Homologous Sequence db
#############################

## extract all the protein sequences into single proteome fasta files
mkdir -p $seqdb
faacomplete=$seqdb/all_complete_proteomes.faa
rm -f $faacomplete*
sed -e 's#\(.\+\)/\([^/]\+$\)#\1/\2/\2_protein\.faa\.gz#g' ${allcomplete}_list > ${faacomplete}_list
for faa in `cat ${faacomplete}_list` ; do zcat $faa >> ${faacomplete} ; echo $faa ; done
grep -c '>' ${faacomplete}
# dereplicate proteins in db
python $dbscripts/dereplicate_fasta.py $faacomplete
echo "$(dateprompt)-- $(grep -c '>' $nrfaacomplete) non-redundant proteins in dataset"

## build database of species-to-sequence
python $dbscripts/allgenome_gff2db.py ${allcomplete}_list $ncbitax/scientific_names.dmp

## clustering of proteome db with  MMSeqs2 
# (https://github.com/soedinglab/MMseqs2,  Steinegger M and Soeding J. Sensitive protein sequence searching for the analysis of massive data sets. bioRxiv, doi: 10.1101/079681 (2016))
# compute the memory use of MMSeq2: M = (7 × N × L + 8 × a^k) bytes, N the number of sequences, L their average size, a the size of the alphabet
mmseqslogs=${entlogs}/mmseqs && mkdir -p $mmseqslogs
# create MMseqs2 db
mmseqs createdb $nrfaacomplete $nrfaacomplete.mmseqs-seqdb &> $mmseqslogs/mmseqs-createdb.log
# perform clustering
mmseqstmp=${enttmp}/mmseqs && mkdir -p $mmseqstmp
families=${seqdb}/protein_families && mkdir -p $families
mmseqsclout=${families}/$(basename $nrfaacomplete).mmseqs-clusterdb_default
# perform similarity search and clustering ; uses all CPU cores by default
mmseqs cluster ${nrfaacomplete}.mmseqs-seqdb $mmseqsclout $mmseqstmp &> $mmseqslogs/$(basename $mmseqsclout).log &
# generate indexed fasta file listing all protein families
mmseqs createseqfiledb ${nrfaacomplete}.mmseqs-seqdb $mmseqsclout ${mmseqsclout}_clusters
# generate separate fasta files with family identifiers distinc from representative sequence name
python $dbscripts/split_mmseqs-clustdb_fasta.py ${mmseqsclout}_clusters "${famprefix}P"
echo "$(dateprompt)-- classified into $(ls ${mmseqsclout}_clusters_fasta/ | wc -l) non-redundant protein clusters"
echo "${datepad}-- including artificial cluster ${famprefix}P000000 gathering $(grep -c '>' ${mmseqsclout}_clusters_fasta/${famprefix}P000000.fasta) ORFan nr proteins"
echo "${datepad}-- (NB: some are not true ORFans as can be be present as identical sequences in several genomes)"

####################################
## 02. Homologous Sequence Alignemnt
####################################

## prepare protein families for alignment
protfamseqs=${mmseqsclout}_clusters_fasta
protali=$entdb/02.clustalo_alignments
nrprotali=$protali/nr_protfam_clustalo_alignments
mkdir -p $nrprotali
python $dbscripts/schedule_ali_task.py $protfamseqs.tab $protfamseqs $protali/$(basename ${protali})_tasklist `nproc`

## align non-redundant protein families
for tasklist in `ls $protali/${protali##*.}_tasklist.*` ; do
  bntl=`basename $tasklist`
  $dbscripts/run_clustalo_sequential.sh $tasklist $nrprotali &> $entlogs/clustalo/$bntl.clustalo.log &
done

## check consistency of non-redundant protein sets
mkdir $raptmp
cd $raptmp
cut -f1 ${allcomplete}_allproteins_info.tab | grep -v "^$\|protein_id" | sort -u > ${allcomplete}_allproteins_in_gff
grep '>' $nrfaacomplete | cut -d' ' -f1 | cut -d'>' -f2 | sort -u > ${nrfaacomplete}_protlist
# compare original dataset of nr protein (as described in the input GFF files) to the aligned nr proteome
diff ${allcomplete}_allproteins_in_gff ${nrfaacomplete}_protlist > $raptmp/diff_prot_info_fasta
if [ -s $raptmp/diff_prot_info_fasta ] ; then 
  >&2 echo "ERROR $(dateprompt): inconsistent propagation of the protein dataset:"
  >&2 echo "present in aligned fasta proteome / absent in info table generated from input GFF:"
  >&2 grep '>' diff_prot_info_fasta | cut -d' ' -f2
  >&2 echo "present in info table generated from input GFF / absent in aligned fasta proteome:"
  >&2 grep '<' diff_prot_info_fasta | cut -d' ' -f2
  exit 1
fi

## reconstruct full (redundant) protein alignments
# make list of CDS sets
for cg in `cat ${allcomplete}_list` ; do ls $cg/*_cds_from_genomic.fna.gz >> ${complete}/all_cds_fasta_list ; done

# generate (full protein alignment, unaligned CDS fasta) file pairs and reverse-translate alignments to get CDS (gene family) alignments
python $dbscripts/extract_full_prot+cds_family_alignments.py ${nrprotali} ${protfamseqs}/${protorfanclust}.fasta ${allcomplete}_allproteins_info.tab \
${allcomplete}_allreplicons_info.tab ${complete}/all_cds_fasta_list $raplogs

#~ ## check consistency of full reverse translated alignement set
ok=1
for fam in `ls /enterobac/PanteroDB_v0.3/02.clustalo_alignments/full_cdsfam_alignments/ | cut -d'.' -f1` ; do 
nseqin1=`grep -c '>' $protali/full_cdsfam_fasta/$fam.fasta`
nseqin2=`grep -c '>' $protali/full_protfam_alignments/$fam.aln`
nseqout=`grep -c '>' $protali/full_cdsfam_alignments/$fam.aln`
if [[ "$nseqin1" != "$nseqout" || "$nseqin2" != "$nseqout" ]] ; then 
  >&2 echo "$fam - full_cdsfam_fasta: ${nseqin1}; full_protfam_alignments: ${nseqin2}; full_cdsfam_alignments: $nseqout"
  ok=0
  echo $fam
fi
done > $raptmp/pal2nal_missed_fams
if [ ok -lt 1 ] ; then
  >&2 echo "WARNING $(dateprompt): failure of pal2nal.pl reverse translation step for families: $(cat $enttmp/pal2nal_missed_fams | xargs)"
  >&2 echo "will use tranposeAlignmentProt2CDS.py instead, a less safe, but more permissive method for generating CDS alignment"
  for fam in `cat $enttmp/pal2nal_missed_fams` ; do
    $dbscripts/tranposeAlignmentProt2CDS.py $protali/full_cdsfam_fasta/$fam.fasta $protali/full_protfam_alignments/$fam.aln $protali/full_cdsfam_alignments/$fam.aln
  done
fi
# protein alignments do not match the CDS sequences
# transpose the CDS into the positions of the aligned protein; assumes no indels, only mismatches and possibly shortenned sequences

# join non-ORFan and ORFan family count matrices
rm -f ${protali}/all_families_genome_counts.mat*
cat ${protali}/full_families_genome_counts-noORFans.mat > ${protali}/all_families_genome_counts.mat
tail -n +2 ${protali}/ENTCGC000000_genome_counts-ORFans.mat >> ${protali}/all_families_genome_counts.mat
gzip ${protali}/all_families_genome_counts.mat


##############################################
## 03. Create and Initiate PostgreSQL database
##############################################

export database=$entdb/03.database
mkdir -p ${database}
cd ${database}
### create PostgreSQL database "panterodb" (pan-genome enterobacteraceae database)
# NB: expect INTERACTTIVE PROMPT here !!
# lower-case name for SQL db: 'panterodb_v0.3'
$dbscripts/create_pangenome_sql_db.sh ${entdbname,,}
