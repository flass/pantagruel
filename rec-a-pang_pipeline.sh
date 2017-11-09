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

#### Set mandatory environment variables / parameters
export myemail='me.myself@respectable-institu.ti.on'
export raproot='/path/to/base/folder/for/all/this/business'
export dbscripts='/path/to/pipeline/scripts'

export rapdbname='aBoldDatabaseName' # mostly name of the top folder
export famprefix='ABCDE'             # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.

## Set facultative parameters (can also be set interactively during the pipeline run)
export pseudocoremingenomes=0		
#~ echo "Enter the minimum number of genomes to possess an unicopy gene family for it to be part of pseudo-core genome:"
#~ while [[ -z $pseudocoremingenomes || $pseudocoremingenomes -lt 1 ]] ; do
 #~ read -p 'Please enter non-null integer value for minimum of genomes represented in pseudo-core unicopy gene families: ' pseudocoremingenomes
#~ done
#~ export pseudocoremingenomes=$pseudocoremingenomes

## load generic environment variables derived from the above
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
ls ${assemblies}/GCF_* -d > ${allcomplete}_withallstrains_list
# complete genomes do not have their sequencing / assembly details reported in *_assembly_stats.txt files ; rather lok for the sequencing metadata in the GenBank flat file
for assemb in `cat ${allcomplete}_withallstrains_list` ; do 
ass=$(basename $assemb) ; tech=`zcat $assemb/${ass}_genomic.gbff.gz | grep -m 1 "Sequencing Technology" | awk -F' :: ' '{print $NF}' | sed -e 's/Il;lumina/Illumina/g'`
if [ -z "$tech" ] ; then tech='NA' ; fi
echo -e "${ass}\t${tech}" ; done > ${indata}/sequencing_technology.tab

## extract assembly/sample metadata from flat files
python $dbscripts/extract_metadata_from_gbff.py --assembly_folder_list=${allcomplete}_withallstrains_list --add_raw_metadata=${complete}/manual_metadata_dictionary.tab \
--add_curated_metadata=${complete}/manual_curated_metadata_dictionary.tab --add_dbxref=${complete}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism"
# remove lab strains
grep "laboratory" ${allcomplete}_withallstrains_metadata_curated.tab | cut -f1 > $complete/labstrains_list
grep "endosymbiont" ${allcomplete}_withallstrains_metadata_curated.tab | cut -f1 > $complete/endosymbiontstrains_list
for strain in `cat $complete/labstrains_list $complete/endosymbiontstrains_list` ; do
 grep $strain ${allcomplete}_withallstrains_list >> ${allcomplete}_excludestrains_list
done
diff ${allcomplete}_withallstrains_list ${allcomplete}_excludestrains_list | grep '<' | cut -d' ' -f2 > ${allcomplete}_list
# regenerate metadata tables without the lab strains
python $dbscripts/extract_metadata_from_gbff.py --assembly_folder_list=${allcomplete}_list --add_raw_metadata=${complete}/manual_metadata_dictionary.tab \
--add_curated_metadata=${complete}/manual_curated_metadata_dictionary.tab --add_dbxref=${complete}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism"

export ngenomes=$((`wc -l ${complete}/complete_genomes_metadata.tab | cut -d' ' -f1` - 1))
echo "work with a database of $ngenomes genomes (excluding lab strains)"

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
python $dbscripts/split_mmseqs-clustdb_fasta.py ${mmseqsclout}_clusters "${famprefix}P" 6
echo "$(dateprompt)-- $(wc -l ${mmseqsclout}_clusters_fasta.tab | cut -d' ' -f1) non-redundant proteins"
echo "$(dateprompt)-- classified into $(ls ${mmseqsclout}_clusters_fasta/ | wc -l) clusters"
echo "${datepad}-- including artificial cluster ${famprefix}P000000 gathering $(grep -c '>' ${mmseqsclout}_clusters_fasta/${famprefix}P000000.fasta) ORFan nr proteins"
echo "${datepad}-- (NB: some are not true ORFans as can be be present as identical sequences in several genomes)"

####################################
## 02. Homologous Sequence Alignemnt
####################################

## prepare protein families for alignment
protfamseqs=${mmseqsclout}_clusters_fasta
export protali=$entdb/02.gene_alignments
nrprotali=$protali/nr_protfam_clustalo_alignments
mkdir -p $nrprotali
# generate lists of many jobs to be executed by a single process (clustalo jobs are too short to be dispatched via a cluster array job)
# clustalomega is already parallelized (for part of the computation, i.e. distance matrix calculation) 
# so less threads than avalilable cores (e.g. use 1 process / 2 cores) are specified
# so each process can run at least at 100% and mobilize 2 cores when parrallel,
# or more when cores unused by other concurent pocesses (during their sequential phase)
python $dbscripts/schedule_ali_task.py $protfamseqs.tab $protfamseqs $protali/$(basename ${protali})_tasklist $((`nproc` / 2))

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
# create database dump
sqldump=${database}/${sqldbname}_dump_$(date +'%d%m%Y')
pg_dump -d ${sqldbname} -j 6 -F d -f $sqldump && chmod -w -R $sqldump/

# generate reference for translation of genome assembly names into short identifier codes (uing UniProt "5-letter" code when available).
psql -d ${sqldbname} -A -t -c "select assembly_id, code from genome.assemblies;" | sed -e 's/ |/\t/g' > $database/genome_codes.tab
# and CDS names into code of type "SPECIES_CDSTAG" (ALE requires such a naming)
psql -d ${sqldbname} -A -t -c "select genbank_cds_id, cds_code from genome.coding_sequences;" | cut -d'|' -f2,3 | sed -e 's/|/\t/g' > $database/cds_codes.tab

# translates the header of the alignment files
export alifastacodedir=${protali}/full_cdsfam_alignments_species_code
mkdir -p ${alifastacodedir}
# use parrallel bash command
$dbscripts/lsfullpath.py ${protali}/full_cdsfam_alignments .aln > ${protali}/full_cdsfam_alignment_list
$dbscripts/genbank2code_fastaseqnames.py ${protali}/full_cdsfam_alignment_list ${database}/cds_codes.tab ${alifastacodedir} > $entdb/logs/genbank2code_fastaseqnames.log


###########################################
## 04. Core Genome Phylogeny (Species tree)
###########################################

# have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
$dbscripts/select_pseudocore_genefams.r 
# new value of pseudocoremingenomes stored in script exit status
export pseudocoremingenomes=$?

## generate core genome alignment path list
export coregenome=$entdb/04.core_genome
mkdir -p ${coregenome}
export pseudocore=${coregenome}/pseudo-core-${pseudocoremingenomes}-unicopy
rm -f ${pseudocore}_cds_aln_list
for fam in `cat ${protali}/pseudo-core-${pseudocoremingenomes}-unicopy_families.tab` ; do ls ${alifastacodedir}/$fam.codes.aln >> ${pseudocore}_cds_aln_list ; done
# concatenate pseudo-core CDS alignments
python $dbscripts/concat.py ${pseudocore}_cds_aln_list ${pseudocore}_concat_cds.aln
rm ${alifastacodedir}/*_all_sp_new

### compute species tree using RAxML
export coretree=${coregenome}/raxml_tree
export treename=${pseudocore}_concat_cds_${ngenomes}entero
mkdir -p ${entlogs}/raxml ${coretree}
# define RAxML options; uses options -T (threads) -j (checkpoints) 
raxmloptions="-n ${treename} -m GTRCATI -j -p 1753 -T 8 -w ${coretree}"
## first check the alignment for duplicate sequences and write a reduced alignment with  will be excluded
raxmlHPC-PTHREADS-AVX -s ${pseudocorealn} ${raxmloptions} -f c &> ${entlogs}/raxml/${treename}.check.log
# 117/880 exactly identical excluded
grep 'exactly identical$' ${coretree}/RAxML_info.${treename} | sed -e 's/IMPORTANT WARNING: Sequences \(.\+\) and \(.\+\) are exactly identical/\1\t\2/g' > ${pseudocorealn}.identical_sequences
rm ${coretree}/RAxML_info.${treename}
## first just single ML tree on reduced alignment
raxmlHPC-PTHREADS-AVX -s ${pseudocorealn}.reduced ${raxmloptions} &> ${entlogs}/raxml/${treename}.ML_run.log &
#~ # resume analysis after crash
#~ raxmlHPC-PTHREADS-AVX -s ${pseudocorealn}.reduced ${raxmloptions} -t ${coretree}/RAxML_checkpoint.${treename}.14 &> ${entlogs}/raxml/${treename}.ML_run.log2 &

## root tree with MAD (Tria et al. Nat Ecol Evol (2017) Phylogenetic rooting using minimal ancestor deviation. s41559-017-0193 doi:10.1038/s41559-017-0193)
R BATCH --vanilla --slave << EOF
source('${dbscripts}/mad_R/MAD.R')
mad.rooting = MAD(readLines('${coretree}/RAxML_bestTree.${treename}'), 'full')
write(mad.rooting[[1]], file='${coretree}/RAxML_bestTree.${treename}.MADrooted')
save(mad.rooting, file='${coretree}/RAxML_bestTree.${treename}.MADrooting.RData')
EOF

## RAxML, with parametric bootstrap
raxmlHPC-PTHREADS-AVX -s ${pseudocorealn}.reduced ${raxmloptions} -b 198237 -N 500 &> ${entlogs}/raxml/${treename}_bs.log &


#############################################################
## 05. Gene trees (full [ML] and collapsed [bayesian sample])
#############################################################

export genetrees=$entdb/05.gene_trees
export mlgenetrees=${genetrees}/raxml_trees
mkdir -p ${mlgenetrees}
mkdir -p $entlogs/raxml/gene_trees

basequery="select gene_family_id, size from genome.gene_family_sizes where gene_family_id is not null and gene_family_id!='$cdsorfanclust'"
python $dbscripts/query_gene_fam_sets.py --db.con.par="dbname=panterodb_v0.3" --outprefix='cdsfams_' --dirout=${protali} \
 --base.query=${basequery} --famsets.min.sizes="4,500,2000,10000" --famsets.max.sizes="499,1999,9999,"

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
allcdsfam2phylo=($(ls ${protali}/cdsfams_*))
allmems=(4 8 32 64)
allwalltimes=(24 24 72 72)
for i in ${!allcdsfam2phylo[@]} ; do
cdsfam2phylo=${allcdsfam2phylo[$i]} ; mem=${allmems[$i]} ; wt=${allwalltimes[$i]}
echo "cdsfam2phylo=${cdsfam2phylo} ; mem_per_core=${mem}gb ; walltime=${wt}:00:00"
tasklist=${cdsfam2phylo}_aln_list
rm -f $tasklist ; for fam in `cut -f1 $cdsfam2phylo` ; do ls ${alifastacodedir}/${fam}.codes.aln >> $tasklist ; done
if [ "$(wc -l $cdsfam2phylo | cut -f1 -d' ')" -lt "$(wc -l $tasklist | cut -f1 -d' ')" ] ; then 
  >&2 echo "ERROR $(dateprompt): missing gene family alignments; please fix the list '$tasklist' or the content of folder '$alifastacodedir/' before continuing."
  exit 1
fi
qsubvars="tasklist=$tasklist,outputdir=$mlgenetrees"
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
chunksize=3000
jobranges=($($dbscripts/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
echo "jobrange=$jobrange"
qsub -J $jobrange -l walltime=${wt}:00:00 -l select=1:ncpus=4:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) -o $entlogs/raxml/gene_trees -j oe -v "$qsubvars" $HOME/scripts/misc/raxml_array_PBS.qsub
done
done

## detect clades to be collapsed in gene trees
export colalinexuscodedir=${protali}/collapsed_cdsfam_alignments_species_code
cladesupp=70
subcladesupp=35
criterion='bs'
withinfun='median'
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
mkdir -p ${colalinexuscodedir}/${collapsecond}
mlgenetreelist=${mlgenetrees%*s}_list
$dbscripts/lsfullpath.py ${mlgenetrees}/${mainresulttag} | sort > ${mlgenetreelist}

Njob=`wc -l ${mlgenetreelist} | cut -f1 -d' '`
chunksize=3000
jobranges=($($dbscripts/get_jobranges.py $chunksize $Njob))

for jobrange in ${jobranges[@]} ; do
beg=`echo $jobrange | cut -d'-' -f1`
tail -n +${beg} ${mlgenetreelist} | head -n ${chunksize} > ${mlgenetreelist}_${jobrange}
qsub -N mark_unresolved_clades -l select=1:ncpus=12:mem=16gb,walltime=24:00:00 -o ${entlogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe -V -S /usr/bin/bash << EOF
module load python
python $dbscripts/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist}_${jobrange} --diraln=${alifastacodedir} --fmt_aln_in='fasta' \
 --threads=12 --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output \
 --clade_stem_conds="[('$criterion', '>=', $cladesupp)]" --within_clade_conds="[('$withinfun', '$criterion', '<=', $subcladesupp, -1), ('max', '$criterion', '<', $cladesupp, -1)]"
EOF
done

## run mrbayes on collapsed alignments
export collapsed_genetrees=${genetrees}/collapsed_mrbayes_trees
mkdir -p ${collapsed_genetrees}/${collapsecond}
alreadytrees=${collapsed_genetrees}_${collapsecond}_list
$dbscripts/lsfullpath.py ${collapsed_genetrees}/${collapsecond} con.tre > $alreadytrees
alreadytasklist=${colalinexuscodedir}/${collapsecond}_ali_list_done
sed -e "s#${collapsed_genetrees}/${collapsecond}/\(.\+\)\.mb\.con\.tre#${colalinexuscodedir}/${collapsecond}/collapsed_alns/\1\.nex#g" $alreadytrees > $alreadytasklist

tasklist=${colalinexuscodedir}/${collapsecond}_ali_list
rm -f $tasklist
$dbscripts/lsfullpath.py ${colalinexuscodedir}/${collapsecond}/collapsed_alns > $tasklist
sort $tasklist > $tasklist.sort
sort $alreadytasklist > $alreadytasklist.sort
dtag=$(date +"%Y-%m-%d-%H-%M-%S")
comm -2 -3  $tasklist.sort $alreadytasklist.sort > ${tasklist}_todo_${dtag}
Njob=`wc -l ${tasklist}_todo | cut -f1 -d' '`
chunksize=1000
jobranges=($($dbscripts/get_jobranges.py $chunksize $Njob))
qsubvar="mbversion=3.2.6, tasklist=${tasklist}_todo_${dtag}, outputdir=${collapsed_genetrees}/${collapsecond}, mbmcmcpopt='Nruns=2 Ngen=2000000 Nchains=6'"
for jobrange in ${jobranges[@]}	; do
echo $jobrange $qsubvar
qsub -J $jobrange -N mb_panterodb -o ${entdlogs}/mrbayes/collapsed_mrbayes_trees_${jobrange}_${dtag} -v "$qsubvar" ${dbscripts}/mrbayes_array_PBS.qsub
done
