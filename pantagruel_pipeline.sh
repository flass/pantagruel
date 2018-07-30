#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 22 August 2017

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

#### Set mandatory environment variables / parameters
export myemail="me.myself@respectable-institu.ti.on"
export raproot="/path/to/base/folder/for/all/this/business"
export ptgrepo="/path/to/pantagruel_repository"
export famprefix="REPLACEfamprefix"
export rapdbname="aBoldDatabaseName" # mostly name of the top folder
export famprefix="ABCDE"             # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.

# derive other important environmnet variables
export ptgscripts="${ptgrepo}/scripts"
export PYTHONPATH=$PYTHONPATH:${ptgscripts}
cd ${ptgrepo} ; export ptgversion=$(git log | grep commit) ; cd -

# create head folders
export rapdb=${raproot}/${rapdbname}
export raplogs=${rapdb}/logs
export raptmp=${rapdb}/tmp
mkdir -p ${rapdb}/ ${raplogs}/ ${raptmp}/

#### Set facultative environment variables / parameters
export pseudocoremingenomes=''       # defaults to empty variable in which case will be set INTERACTIVELY at stage 04.core_genome of the pipeline
envsourcescript=${HOME}/environ_pantagruel_${rapdbname}.sh

rm -f ${raptmp}/sedenvvar.sh
echo -n "cat ${ptgscripts}/environ_pantagruel_template.sh" > ${raptmp}/sedenvvar.sh
for var in myemail raproot ptgscripts famprefix rapdbname famprefix pseudocoremingenomes ; do
echo -n " | sed -e \"s#REPLACE${var}#${var}#\"" >> ${raptmp}/sedenvvar.sh
done
echo -n " > ${envsourcescript}" >> ${raptmp}/sedenvvar.sh
bash < ${raptmp}/sedenvvar.sh
## load generic environment variables derived from the above
source ${envsourcescript}
cat "source ${envsourcescript}" >> ~/.bashrc

# folders for optional custom genomes
export customassemb=${raproot}/user_genomes
mkdir -p ${customassemb}/contigs/
echo "sequencing.project.id,genus,species,strain,taxid,locus_tag_prefix" | tr ',' '\t' > ${straininfo}

echo "please copy/link raw sequence (in multi-fasta format) files of custom (user-generated) assemblies into ${customassemb}/contigs/"
echo "and fill up the information table ${straininfo} (tab-delimited fields) according to header:"
cat ${straininfo}

#################################
## 00. Data download and grooming
#################################

### Download data
## using FTP
# downloading relevant tables from NCBI Taxonomy database 
user='anonymous'
pswd=$myemail
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
source="/pub/taxonomy"
files="taxcat.tar.gz* taxcat_readme.txt taxdump.tar.gz* taxdump_readme.txt"
ncbitax="${raproot}/NCBI/Taxonomy"
mkdir -p $ncbitax
lftp -c "open $openparam; cd $source ; mget -O $ncbitax/ $files ; quit"
mv $ncbitax/README $ncbitax/readme_accession2taxid
# reduce taxonomic id2name reference file complexity
for tgz in `ls $ncbitax/*.tar.gz` ; do md5sum -c $tgz.md5 && tar -xzf $tgz ; done

## using NCBI web interface:
# Search Assembly database using a query defined as convenient:	
# e.g.: 'txid543[Organism:exp] AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter]) AND ("complete genome"[filter])'
# save assemblies (Download Assemblies > Source Database: RefSeq; File type: All file type (including assembly-structure directory)) as
$indata/genome_assemblies.tar
# extract assembly data
tar -C $indata/assemblies/ -xf genome_assemblies.tar
assfolder=`ls -d $indata/ncbi-genomes-* | head -n 1`
rm genome_assemblies.tar 
# store in centralised folder for NCBI assemblies and just record links
mv $assfolder/* $ncbiass/
mv $assfolder/ $indata/assemblies/
for ass in `cut -f3 $indata/assembly_result.txt | tail -n +2` ; do ln -s $ncbiass/${ass}_* $indata/assemblies/  ; done


### can be automated using NCBI Entrez tools
#~ example from https://www.ncbi.nlm.nih.gov/books/NBK179288/ :
 #~ for org in \
    #~ "Agrobacterium tumefaciens" \
    #~ "Bacillus anthracis" \
    #~ "Escherichia coli" \
    #~ "Neisseria gonorrhoeae" \
    #~ "Pseudomonas aeruginosa" \
    #~ "Shigella flexneri" \
    #~ "Streptococcus pneumoniae"
  #~ do
    #~ esearch -db assembly -query "$org [ORGN]" |
    #~ efilter -query "representative [PROP]" |
    #~ elink -target nuccore -name assembly_nuccore_refseq |
    #~ efetch -format docsum |
    #~ xtract -pattern DocumentSummary -element AccessionVersion Slen Title |
    #~ sed 's/,.*//' |
    #~ grep -v -i -e scaffold -e contig -e plasmid -e sequence -e patch |
    #~ sort -t $'\t' -k 2,2nr
  #~ done

## if the folder of custom/user-provided set of genomes is not empty
if [[ "$(ls -A "${customassemb}/contigs/" 2>/dev/null)" ]] ; then

### uniformly annotate the raw sequence data with Prokka
# first add all the (representative) proteins in the dataset to the custom reference prot database for Prokka to search for similarities
ls ${indata}/assemblies/*/*gbff.gz > ${indata}/assemblies_gbffgz_list
parallel -a ${indata}/assemblies_gbffgz_list 'gunzip -k'
# extract protein sequences
prokka-genbank_to_fasta_db ${indata}/assemblies/*/*gbff > ${raptmp}/${refgenus}.faa >& ${raptmp}/prokka-genbank_to_fasta_db.log
# cluster similar sequences
cdhit -i ${raptmp}/${refgenus}.faa -o ${raptmp}/${refgenus}_representative.faa -T 0 -M 0 -G 1 -s 0.8 -c 0.9 &> cdhit.log
# 805428 representative proteins
rm -fv ${raptmp}/${refgenus}.faa ${raptmp}/${refgenus}_representative.faa.clstr
# replace database name for genus detection by Prokka 
prokkablastdb=$(dirname $(dirname $(ls -l `which prokka` | awk '{ print $NF }')))/db/genus/
cp -p ${raptmp}/${refgenus}_representative.faa ${prokkablastdb}/${refgenus}
cd ${prokkablastdb}/
makeblastdb -dbtype prot -in ${refgenus}
cd -

# run Prokka 
export annot=${customassemb}/prokka_annotation
export straininfo=${customassemb}/strain_infos.txt
mkdir -p $annot
for allcontigs in `ls ${customassemb}/contigs/` ; do
gproject=${allcontigs%_assembly*}
date
echo "### assembly: $gproject; contig files from: ${customassemb}/contigs/${allcontigs}"
echo "running Prokka..."
$ptgscripts/run_prokka.sh ${gproject} ${customassemb}/contigs/${allcontigs} ${straininfo} ${annot}/${gproject} ${seqcentre}
echo "fix prokka output to integrate region information into GFF files"
annotoutgff=$(ls ${annot}/${gproject}/*.gff)
python $ptgscripts/add_region_feature2prokkaGFF.py ${annotoutgff[0]} ${annotoutgff[0]%*.gff}.ptg.gff ${straininfo} ${customassemb}/contigs/${allcontigs} 'Unicycler:0.4.3'
echo "fix prokka output to integrate taxid information into GBK files"
annotoutgbk=$(ls ${annot}/${gproject}/*.gbk)
python $ptgscripts/add_taxid_feature2prokkaGBK.py ${annotoutgbk[0]} ${annotoutgbk[0]%*.gbk}.ptg.gbk ${straininfo}
echo "done."
date
done &> $raptmp/prokka.log &

# create assembly directory and link/create relevant files
export gblikeass=${customassemb}/genbank-format_assemblies
mkdir -p $gblikeass/
for gproject in `ls $annot` ; do
gff=$(ls ${annot}/${gproject}/ | grep 'ptg.gff')
assemb=${gproject}.1_${gff[0]%*.ptg.gff}
assembpathdir=${gblikeass}/${assemb}
assembpathrad=${assembpathdir}/${assemb}
mkdir -p ${assembpathdir}/
ls ${assembpathdir}/ -d
gffgz=${assembpathrad}_genomic.gff.gz
gzip -c ${annot}/${gproject}/$gff > $gffgz
ls $gffgz
gbk=$(ls ${annot}/${gproject}/ | grep 'ptg.gbk') ; gbkgz=${assembpathrad}_genomic.gbff.gz
gzip -c ${annot}/${gproject}/$gbk > $gbkgz
ls $gbkgz
faa=$(ls ${annot}/${gproject}/ | grep '.faa') ; faagz=${assembpathrad}_protein.faa.gz
gzip -c ${annot}/${gproject}/$faa > $faagz
ffn=$(ls ${annot}/${gproject}/ | grep '.ffn') ; fnagz=${assembpathrad}_cds_from_genomic.fna.gz
python ${ptgscripts}/format_prokkaCDS.py ${annot}/${gproject}/$ffn $fnagz
ls $fnagz
done

ln -s ${gblikeass}/* ${indata}/assemblies/
fi
#### end of the block treating custom genome set

### Groom data
## generate assembly statistics to verify genome finishing status
${ptgscripts}/groom_refseq_data.sh ${indata}/assemblies ${indata}/assembly_stats

manuin=${genomeinfo}/manual_input_metadata
mkdir ${genomeinfo}/metadata_${rapdbname}/ ${manuin}/
touch ${manuin}/manual_metadata_dictionary.tab ${manuin}/manual_curated_metadata_dictionary.tab ${manuin}/manual_dbxrefs.tab
rm -rf ${genomeinfo}/assemblies*
ls ${assemblies}/* -d > ${genomeinfo}/assemblies_list
mkdir -p ${genomeinfo}/assembly_metadata
## extract assembly/sample metadata from flat files
python ${ptgscripts}/extract_metadata_from_gbff.py --assembly_folder_list=${genomeinfo}/assemblies_list --add_raw_metadata=${manuin}/manual_metadata_dictionary.tab \
--add_curated_metadata=${manuin}/manual_curated_metadata_dictionary.tab --add_dbxref=${manuin}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism" --output=${genomeinfo}/assembly_metadata

export ngenomes=$((`wc -l ${genomeinfo}/assembly_metadata/metadata.tab | cut -d' ' -f1` - 1))
echo "work with a database of $ngenomes genomes (excluding lab strains)"

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
mmseqslogs=${raplogs}/mmseqs && mkdir -p $mmseqslogs
# create MMseqs2 db
mmseqs createdb ${allfaarad}.nr.faa ${allfaarad}.nr.mmseqsdb &> $mmseqslogs/mmseqs-createdb.log
# perform clustering
mmseqstmp=${raptmp}/mmseqs && rm -rf $mmseqstmp && mkdir -p $mmseqstmp
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

####################################
## 02. Homologous Sequence Alignemnt
####################################

## prepare protein families for alignment
export protfamseqs=${mmseqsclout}_clusters_fasta
export protali=${rapdb}/02.gene_alignments
nrprotali=$protali/nr_protfam_clustalo_alignments
mkdir -p ${nrprotali}/ ${raplogs}/clustalo/
# generate lists of many jobs to be executed by a single process (clustalo jobs are too short to be dispatched via a cluster array job)
# clustalomega is already parallelized (for part of the computation, i.e. distance matrix calculation) 
# so less threads than avalilable cores (e.g. use 1 process / 2 cores) are specified
# so each process can run at least at 100% and mobilize 2 cores when parrallel,
# or more when cores unused by other concurent pocesses (during their sequential phase)
python ${ptgscripts}/schedule_ali_task.py $protfamseqs.tab $protfamseqs $protali/$(basename ${protali})_tasklist $((`nproc` / 2))

## align non-redundant protein families
for tasklist in `ls $protali/${protali##*.}_tasklist.*` ; do
  bntl=`basename $tasklist`
  ${ptgscripts}/run_clustalo_sequential.sh $tasklist $nrprotali 2 &> $raplogs/clustalo/$bntl.clustalo.log &
done

## check consistency of non-redundant protein sets
mkdir -p $raptmp
protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'nr_protein_id' | cut -d':' -f1)
if [ -z $protidfield ] ; then 
 protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'protein_id' | cut -d':' -f1)
fi
cut -f $protidfield ${genomeinfo}/assembly_info/allproteins_info.tab | grep -v "^$\|protein_id" | sort -u > ${genomeinfo}/assembly_info/allproteins_in_gff
grep '>' ${allfaarad}.nr.faa | cut -d' ' -f1 | cut -d'>' -f2 | sort -u > ${allfaarad}.nr_protlist
# compare original dataset of nr protein (as described in the input GFF files) to the aligned nr proteome
diff ${genomeinfo}/assembly_info/allproteins_in_gff ${allfaarad}.nr_protlist > $raptmp/diff_prot_info_fasta
if [ -s $raptmp/diff_prot_info_fasta ] ; then 
  >&2 echo "ERROR $(dateprompt): inconsistent propagation of the protein dataset:"
  >&2 echo "present in aligned fasta proteome / absent in info table generated from input GFF:"
  >&2 grep '>' $raptmp/diff_prot_info_fasta | cut -d' ' -f2
  >&2 echo "present in info table generated from input GFF / absent in aligned fasta proteome:"
  >&2 grep '<' $raptmp/diff_prot_info_fasta | cut -d' ' -f2
  exit 1
fi

## reconstruct full (redundant) protein alignments
# make list of CDS sets
for cg in `cat ${genomeinfo}/assemblies_list` ; do ls $cg/*_cds_from_genomic.fna.gz >> ${genomeinfo}/all_cds_fasta_list ; done

# generate (full protein alignment, unaligned CDS fasta) file pairs and reverse-translate alignments to get CDS (gene family) alignments
mkdir -p ${raplogs}/extract_full_prot_and_cds_family_alignments/
python ${ptgscripts}/extract_full_prot_and_cds_family_alignments.py --nrprot_fam_alns ${nrprotali} --singletons ${protfamseqs}/${protorfanclust}.fasta \
 --prot_info ${genomeinfo}/assembly_info/allproteins_info.tab --repli_info ${genomeinfo}/assembly_info/allreplicons_info.tab --assemblies ${assemblies} \
 --dirout ${protali} --famprefix ${famprefix} --logs ${raplogs}/extract_full_prot_and_cds_family_alignments --identical_prots ${allfaarad}.identicals.list

## check consistency of full reverse translated alignment set
ok=1
for fam in `ls $protali/full_cdsfam_alignments/ | cut -d'.' -f1` ; do 
nseqin1=`grep -c '>' $protali/full_cdsfam_fasta/$fam.fasta`
nseqin2=`grep -c '>' $protali/full_protfam_alignments/$fam.aln`
nseqout=`grep -c '>' $protali/full_cdsfam_alignments/$fam.aln`
if [[ "$nseqin1" != "$nseqout" || "$nseqin2" != "$nseqout" ]] ; then 
  >&2 echo "$fam - full_cdsfam_fasta: ${nseqin1}; full_protfam_alignments: ${nseqin2}; full_cdsfam_alignments: $nseqout"
  ok=0
  echo $fam
fi
done > $raptmp/pal2nal_missed_fams
if [ $ok -lt 1 ] ; then
  >&2 echo "WARNING $(dateprompt): failure of pal2nal.pl reverse translation step for families: $(cat $raptmp/pal2nal_missed_fams | xargs)"
  >&2 echo "will use tranposeAlignmentProt2CDS.py instead, a less safe, but more permissive method for generating CDS alignment"
  for fam in `cat $raptmp/pal2nal_missed_fams` ; do
    ${ptgscripts}/tranposeAlignmentProt2CDS.py $protali/full_cdsfam_fasta/$fam.fasta $protali/full_protfam_alignments/$fam.aln $protali/full_cdsfam_alignments/$fam.aln
  done
fi
# protein alignments do not match the CDS sequences
# transpose the CDS into the positions of the aligned protein; assumes no indels, only mismatches and possibly shortenned sequences

# join non-ORFan and ORFan family count matrices
rm -f ${protali}/all_families_genome_counts.mat*
cat ${protali}/full_families_genome_counts-noORFans.mat > ${protali}/all_families_genome_counts.mat
tail -n +2 ${protali}/${famprefix}C000000_genome_counts-ORFans.mat >> ${protali}/all_families_genome_counts.mat
gzip ${protali}/all_families_genome_counts.mat

##############################################
## 03. Create and Populate SQLite database
##############################################

export database=${rapdb}/03.database
export sqldbname=${rapdbname,,}
export sqldb=${database}/${sqldbname}
mkdir -p ${database}
cd ${database}
### create and populate SQLite database

## Genome schema: initiate and populate
${ptgscripts}/pantagruel_sqlite_genome_db.sh ${database} ${sqldbname} ${genomeinfo}/metadata_${rapdbname} ${genomeinfo}/assembly_info ${protali} ${protfamseqs}.tab ${protorfanclust} ${cdsorfanclust} ${straininfo}

# dump reference table for translation of genome assembly names into short identifier codes (uing UniProt "5-letter" code when available).
sqlite3 ${sqldb} "select assembly_id, code from assemblies;" | sed -e 's/|/\t/g' > ${database}/genome_codes.tab
sqlite3 ${sqldb} "select code, organism, strain from assemblies;" | sed -e 's/|/\t/g' > ${database}/organism_codes.tab
# and for CDS names into code of type "SPECIES_CDSTAG" (ALE requires such a naming)
# here split the lines with '|' and remove GenBank CDS prefix that is always 'lcl|'
sqlite3 ${sqldb} "select genbank_cds_id, cds_code from coding_sequences;" | sed -e 's/lcl|//g'  | sed -e 's/|/\t/g' > ${database}/cds_codes.tab

# translates the header of the alignment files
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
for mol in prot cds ; do
  alifastacodedir=${protali}/full_${mol}fam_alignments_species_code
  mkdir -p ${alifastacodedir}
  eval "export ${mol}alifastacodedir=${alifastacodedir}"
  # use multiprocessing python script
  ${ptgscripts}/lsfullpath.py ${protali}/full_${mol}fam_alignments .aln > ${protali}/full_${mol}fam_alignment_list
  ${ptgscripts}/genbank2code_fastaseqnames.py ${protali}/full_${mol}fam_alignment_list ${database}/cds_codes.tab ${alifastacodedir} > ${rapdb}/logs/genbank2code_fastaseqnames.${mol}.log
done

## Phylogeny schema: initiate
sqlite3 ${dbfile} < $ptgscripts/pantagruel_sqlitedb_phylogeny_initiate.sql

###########################################
## 04. Core Genome Phylogeny (Species tree)
###########################################

## select psedo-core genome marker gene set
export coregenome=${rapdb}/04.core_genome
mkdir -p ${coregenome}/
#~ if [ -z ${pseudocoremingenomes} ] ; then
  #~ mkdir -p ${coregenome}/pseudo-coregenome_sets/
  #~ # have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
  #~ ${ptgscripts}/select_pseudocore_genefams.r ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets
  #~ # new value of pseudocoremingenomes stored in script exit status
  #~ export pseudocoremingenomes=$?
#~ fi
unset pseudocoremingenomes
mkdir -p ${coregenome}/pseudo-coregenome_sets/
# have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
let "t = ($ngenomes * 3 / 4)" ; let "u = $t - ($t%20)" ; seq $u 20 $ngenomes | cat > ${raptmp}/mingenom ; echo "0" >> ${raptmp}/mingenom
Rscript --vanilla --silent ${ptgscripts}/select_pseudocore_genefams.r ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets < ${raptmp}/mingenom

mv ${rapdb}/environ_pantagruel_${rapdbname}.sh ${rapdb}/environ_pantagruel_${rapdbname}.sh0 && \
 sed -e "s#'REPLACEpseudocoremingenomes'#$pseudocoremingenomes#" ${rapdb}/environ_pantagruel_${rapdbname}.sh0 > ${rapdb}/environ_pantagruel_${rapdbname}.sh
 
## generate core genome alignment path list
export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
rm -f ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list
for fam in `cat ${coregenome}/pseudo-coregenome_sets/${pseudocore}_families.tab` ; do
 ls ${protalifastacodedir}/$fam.codes.aln >> ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list
done
export pseudocorealn=${coregenome}/${pseudocore}_concat_prot.aln
# concatenate pseudo-core prot alignments
python ${ptgscripts}/concat.py ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list ${pseudocorealn}
rm ${protalifastacodedir}/*_all_sp_new

### compute species tree using RAxML
export coretree=${coregenome}/raxml_tree
export treename=${pseudocore}_concat_prot_${ngenomes}-genomes_${rapdbname}
mkdir -p ${raplogs}/raxml ${coretree}
# define RAxML binary and options; uses options -T (threads) -j (checkpoints) 
if [ ! -z $(grep -o avx2 /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-AVX2 -U'
elif [ ! -z $(grep -o avx /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-AVX -U'
elif [ ! -z $(grep -o sse3 /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-SSE3 -U'
else
  raxmlflav=''
fi
raxmlbin="raxmlHPC-PTHREADS${raxmlflav} -T $(nproc)"
raxmloptions="-n ${treename} -m PROTCATLGX -j -p 1753 -w ${coretree}"

## first check the alignment for duplicate sequences and write a reduced alignment with  will be excluded
$raxmlbin -s ${pseudocorealn} ${raxmloptions} -f c &> ${raplogs}/raxml/${treename}.check.log
# 117/880 exactly identical excluded
grep 'exactly identical$' ${coretree}/RAxML_info.${treename} | sed -e 's/IMPORTANT WARNING: Sequences \(.\+\) and \(.\+\) are exactly identical/\1\t\2/g' > ${pseudocorealn}.identical_sequences
rm ${coretree}/RAxML_info.${treename}
## first just single ML tree on reduced alignment
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} &> ${raplogs}/raxml/${treename}.ML_run.log
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bestTree.${treename}

## RAxML, with rapid bootstrap
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -x 198237 -N 1000 &> ${raplogs}/raxml/${treename}_bs.log
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bootstrap.${treename}
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -f b -z ${coretree}/RAxML_bootstrap.${treename} -t ${coretree}/RAxML_bestTree.${treename} 
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bipartitions.${treename}
#~ ## root ML tree
#~ $raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -f I -t RAxML_bipartitionsBranchLabels.${treename}  
#~ mv RAxML_info.${treename} RAxML_info_rootedTree.${treename}

## root tree with MAD (Tria et al. Nat Ecol Evol (2017) Phylogenetic rooting using minimal ancestor deviation. s41559-017-0193 doi:10.1038/s41559-017-0193)
nrbesttree=${coretree}/RAxML_bestTree.${treename}
rootingmethod='MAD'
nrrootedtree=${nrbesttree}.${rootingmethod}rooted
R BATCH --vanilla --slave << EOF
source('${ptgscripts}/mad_R/MAD.R')
mad.rooting = MAD(readLines('${nrbesttree}'), 'full')
write(mad.rooting[[1]], file='${nrrootedtree}')
save(mad.rooting, file='${nrbesttree}.${rootingmethod}rooting.RData}')
EOF

# reintroduce excluded species name into the tree + root as MAD-rooted species tree
nrbiparts=${nrbesttree/bestTree/bipartitions}
export speciestree=${nrrootedtree}.full
python ${ptgscripts}/putidentseqbackintree.py --input.nr.tree=${nrbiparts} --ref.rooted.nr.tree=${nrrootedtree} \
 --list.identical.seq=${pseudocorealn}.identical_sequences --output.tree=${speciestree}

python ${ptgscripts}/code2orgnames_in_tree.py ${speciestree} $database/organism_codes.tab ${speciestree}.names

### generate ultrametric 'dated' species tree for more meaningfulgraphical representation of reconciliation AND to use the dated (original) version of ALE
## use LSD (v 0.3 beta; To et al. Syst. Biol. 2015) to generate an ultrametric tree from the rooted ML tree (only assumin different, uncorrelated clocks for each branch)
alnlen=$( head -n1 ${pseudocorealn}.reduced | cut -d' ' -f2 )
lsd -i ${speciestree} -c -v 1 -s $alnlen -o ${speciestree}.lsd

## delineate populations of near-identical strains (based on tree with branch lengths in subs/site) and generate the population tree, i.e. the species tree withs population collapsed
python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -S ${speciestree} --threads=8 \
--pop_stem_conds="[('lg', '>=', 0.0005), ('bs', '>=', 80)]" --within_pop_conds="[('max', 'lg', '<=', 0.0005, -1)]"
nspepop=$(tail -n+3 ${speciestree%.*}_populations | wc -l)
echo "Found $nspepop disctinct populations in the species tree"


## extract sister clade pairs from the reference species tree for later clade-specific gene search
python << EOF
import tree2
reftree = tree2.Node(file='${speciestree}')
nfout = "${speciestree}_clade_defs"
fout = open(nfout, 'w')
fout.write('\t'.join(['', 'clade', 'sisterclade'])+'\n')
k = 0
for node in reftree:
  if len(node.children) != 2: continue
  for child in node.children:
    if child.nb_leaves() <= 1: continue
    focchildlab = child.label()
    if not focchildlab:
      focchildlab = "clade%d"%k
      k += 1
    focchildleaflabset = ','.join(sorted(child.get_leaf_labels()))
    sischildleaflabset = ','.join(sorted(child.go_brother().get_leaf_labels()))
    fout.write('\t'.join([focchildlab, focchildleaflabset, sischildleaflabset])+'\n')

fout.close()
EOF


#############################################################
## 05. Gene trees (full [ML] and collapsed [bayesian sample])
#############################################################

export genetrees=${rapdb}/05.gene_trees
export mlgenetrees=${genetrees}/raxml_trees
mkdir -p ${mlgenetrees}
mkdir -p $raplogs/raxml/gene_trees

basequery="select gene_family_id, size from gene_family_sizes where gene_family_id is not null and gene_family_id!='$cdsorfanclust'"
python ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${protali} \
 --base.query="${basequery}" --famsets.min.sizes="4,500,2000,10000" --famsets.max.sizes="499,1999,9999,"

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
allcdsfam2phylo=($(ls ${protali}/cdsfams_*))
allmems=(4 8 32 64)
allwalltimes=(24 24 72 72)
allncpus=(4 4 4 16)
for i in ${!allcdsfam2phylo[@]} ; do
cdsfam2phylo=${allcdsfam2phylo[$i]} ; mem=${allmems[$i]} ; wt=${allwalltimes[$i]} ; ncpus=${allncpus[$i]}
echo "cdsfam2phylo=${cdsfam2phylo} ; mem_per_core=${mem}gb ; walltime=${wt}:00:00 ; ncpus=${ncpus}"
tasklist=${cdsfam2phylo}_aln_list
rm -f $tasklist ; for fam in `cut -f1 $cdsfam2phylo` ; do ls ${cdsalifastacodedir}/${fam}.codes.aln >> $tasklist ; done
if [ "$(wc -l $cdsfam2phylo | cut -f1 -d' ')" -lt "$(wc -l $tasklist | cut -f1 -d' ')" ] ; then 
  >&2 echo "ERROR $(dateprompt): missing gene family alignments; please fix the list '$tasklist' or the content of folder '$alifastacodedir/' before continuing."
  exit 1
fi
qsubvars="tasklist=$tasklist,outputdir=$mlgenetrees,reducedaln=true,nbthreads=${ncpus}"
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
# accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
chunksize=3000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
echo "jobrange=$jobrange"
qsub -J $jobrange -l walltime=${wt}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) -o $raplogs/raxml/gene_trees -j oe -v "$qsubvars" ${ptgscripts}/raxml_array_PBS.qsub
done
done


#### OPTION: edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
if [ -z collapseCladeOptions ] ; then
  chaintype='fullgenetree'
  export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
  # not implemented yet
  # only have to only convert the alignments from fasta to nexus
  for aln in `ls ${alifastacodedir}` ; do
    convalign  -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  export nexusaln4chains=${colalinexuscodedir}
  export mboutputdir=${bayesgenetrees}
  
else
  chaintype='collapsed'
  if [-z ${collapsecolid} ] ; then
    collapsecolid=1
  fi
  eval "$collapseCladeOptions"
  # e.g.:  eval "cladesupp=70 ; subcladesupp=35 ; criterion='bs' ; withinfun='median'"

  ## detect clades to be collapsed in gene trees
  export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
  export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
  export collapsecriteriondef="--clade_stem_conds=\"[('$criterion', '>=', $cladesupp)]\" --within_clade_conds=\"[('$withinfun', '$criterion', '<=', $subcladesupp, -1), ('max', '$criterion', '<', $cladesupp, -1)]\""
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  mlgenetreelist=${mlgenetrees%*s}_list
  ${ptgscripts}/lsfullpath.py ${mlgenetrees}/${mainresulttag} | sort > ${mlgenetreelist}
  
  # accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
  Njob=`wc -l ${mlgenetreelist} | cut -f1 -d' '`
  chunksize=3000
  jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
  ncpus=4

  for jobrange in ${jobranges[@]} ; do
  beg=`echo $jobrange | cut -d'-' -f1`
  tail -n +${beg} ${mlgenetreelist} | head -n ${chunksize} > ${mlgenetreelist}_${jobrange}
  qsub -N mark_unresolved_clades -l select=1:ncpus=${ncpus}:mem=16gb,walltime=4:00:00 -o ${raplogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe -V -S /usr/bin/bash << EOF
  module load python
  python ${ptgscripts}/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist}_${jobrange} --diraln=${alifastacodedir} --fmt_aln_in='fasta' \
   --threads=${ncpus} --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output --dir_identseq=${mlgenetrees}/identical_sequences \
   ${collapsecriteriondef}
EOF
  done
  export collapsecoldate=$(date +%Y-%m-%d)
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
  export mboutputdir=${bayesgenetrees}/${collapsecond}
  
fi
#### end OPTION

## run mrbayes on collapsed alignments
export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
mkdir -p ${mboutputdir}
nchains=4
nruns=2
ncpus=$(( $nchains * $nruns ))
tasklist=${nexusaln4chains}_ali_list
rm -f $tasklist
${ptgscripts}/lsfullpath.py ${nexusaln4chains} > $tasklist

#~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
#~ alreadytrees=${mboutputdir}_list
#~ ${ptgscripts}/lsfullpath.py ${mboutputdir} con.tre > $alreadytrees
#~ alreadytasklist=${nexusaln4chains}_ali_list_done
#~ sed -e "s#${mboutputdir}/\(.\+\)\.mb\.con\.tre#${nexusaln4chains}/\1\.nex#g" $alreadytrees > $alreadytasklist
#~ sort $tasklist > $tasklist.sort
#~ sort $alreadytasklist > $alreadytasklist.sort
#~ dtag=$(date +"%Y-%m-%d-%H-%M-%S")
#~ comm -2 -3  $tasklist.sort $alreadytasklist.sort > ${tasklist}_todo_${dtag}
#~ Njob=`wc -l ${tasklist}_todo_${dtag} | cut -f1 -d' '`
#~ qsubvar="mbversion=3.2.6, tasklist=${tasklist}_todo_${dtag}, outputdir=${mboutputdir}, mbmcmcpopt='Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'"
# otherwise could just use $tasklist
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
chunksize=1000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
qsubvar="mbversion=3.2.6, tasklist=${tasklist}, outputdir=${mboutputdir}, mbmcmcpopt='Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'"
for jobrange in ${jobranges[@]} ; do
 echo $jobrange $qsubvar
 qsub -J $jobrange -N mb_panterodb -l select=1:ncpus=${ncpus}:mem=16gb -o ${raplogs}/mrbayes/collapsed_mrbayes_trees_${dtag}_${jobrange} -v "$qsubvar" ${ptgscripts}/mrbayes_array_PBS.qsub
done

#### OPTION: were the rake lades in gene trees collapsed? if yes, these need to be replaced by mock population leaves
if [ -z collapseCladeOptions ] ; then
  export chaintype='fullgenetree'
  export coltreechains=${genetrees}/${chaintype}_tree_chains
  # not implemented yet
  # must generalize the script to only convert the tree chains, not replacing anything in them
  #~ python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -o ${coltreechains} --threads=${ncpus} --reuse=0 --verbose=0 --logfile=${repllogs}_${replrun}.log &

  else
    if [-z ${replacecolid} ] ; then
     replacecolid=1
    fi
  export chaintype='collapsed'
  export coltreechains=${genetrees}/${chaintype}_tree_chains
  export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'
  mkdir -p ${coltreechains}/${collapsecond}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  mkdir -p ${rapdb}/logs/replspebypop
  tasklist=${coltreechains}_${collapsecond}_nexus_list
  ls $bayesgenetrees/${collapsecond}/*run1.t > $tasklist
  repllogd=${rapdb}/logs/replspebypop
  repllogs=$repllogd/replace_species_by_pop_in_gene_trees
  replrun=$(date +'%d%m%Y')

  # local parallel run
  python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.newick -o ${coltreechains}/${collapsecond} \
   --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk --population_node_distance=${speciestree%.*}_interNodeDistPopulations \
   --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${colmethod} --threads=${ncpus} --reuse=0 --verbose=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log &

  ## OR

  #~ # PBS-submitted parallel job
  #~ qsub -N replSpePopinGs -l select=1:ncpus=${ncpus}:mem=64gb,walltime=24:00:00 -o $repllogd -j oe -V << EOF
  #~ module load python
  #~ python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.newick -o ${coltreechains}/${collapsecond} \
   #~ --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk --population_node_distance=${speciestree%.*}_interNodeDistPopulations \
   #~ --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${colmethod} --threads=${ncpus} --reuse=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log
  #~ EOF

  export replacecoldate=$(date +%Y-%m-%d)

  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${dbfile} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} ${collapsecriteriondef} ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}

fi
#### end OPTION


#############################################
## 06. gene tree/ Species tree reconciliations
#############################################

export alerec=${rapdb}/06.ALE_reconciliation
mkdir -p ${alerec}

### perform reconciliations with ALE

# parameters to be set
export ALEversion='v0.4'
export ALEalgo='ALEml_undated'
export recsamplesize=1000
export ALEsourcenote='program compiled from source code from of https://github.com/ssolo/ALE/commits/63f0a3c964074a15f61fd45156ab9e10b5dd45ef'
if [-z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${ALEalgo} == 'ALEml_undated' ] ; then
  export rectype='undat'
else
  export rectype='dated'
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_recs

tasklist=${coltreechains}/${collapsecond}/${colmethod}_Gtrees_list
ls ${coltreechains}/${collapsecond}/${colmethod}/*-Gtrees.nwk > $tasklist
alelogs=${rapdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${colmethod}/${reccol}
mkdir -p $outrecdir

Njob=`wc -l $tasklist | cut -f1 -d' '`
qsubvars="tasklist=$tasklist, resultdir=$outrecdir, spetree=Stree.nwk, nrecs=${recsamplesize}, alealgo=${ALEalgo}"
qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=1:mem=20gb,walltime=24:00:00 -o $alelogs/${reccol} -j oe -v "$qsubvars" ${ptgscripts}/ALE_array_PBS.qsub

export reccoldate=$(date +%Y-%m-%d)

### parse the inferred scenarios
# parameters to be set
export evtypeparse='ST'
export minevfreqparse=0.1
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi
# derived parameters
export parsedreccol=${reccol}_parsed_${parsedreccolid}
export parsedrecs=${alerec}/parsed_recs/${parsedreccol}

mkdir -p ${parsedrecs}
reclist=$outrecdir/ale_collapsed_undat_uml_rec_list
${ptgscripts}/lsfullpath.py $outrecdir/ale_collapsed_undat uml_rec > $reclist
 
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families
python ${ptgscripts}/parse_collapsedALE_scenarios.py --rec_sample_list ${reclist} \
 --populations ${speciestreeBS%.*}_populations --reftree ${speciestreeBS}.lsd.newick \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads 8  &> $entlogs/parse_collapsedALE_scenarios.log &

export parsedreccoldate=$(date +%Y-%m-%d)

## store reconciliation parameters and load parsed reconciliation data into database
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh ${database} ${sqldb} ${parsedrecs} ${ALEversion} ${ALEalgo} ${ALEsourcenote} ${parsedreccol} ${parsedreccolid} ${parsedreccoldate}

# rapid survey of event density over the reference tree
for freqthresh in 0.1 0.25 0.5 ; do
sqlite3 ${sqldb} """
.mode tabs 
select don_branch_id, don_branch_name, rec_branch_id, rec_branch_name, event_type, nb_lineages, cum_freq, cum_freq/nb_lineages as avg_support from (
 select don_branch_id, don_stree.branch_name as don_branch_name, rec_branch_id, rec_stree.branch_name as rec_branch_name, event_type, count(*) as nb_lineages, sum(freq)::real/${nsample} as cum_freq
  from gene_lineage_events 
  inner join species_tree_events using (event_id) 
  inner join species_tree as rec_stree on rec_branch_id=rec_stree.branch_id
  left join species_tree as don_stree on don_branch_id=don_stree.branch_id
 where freq >= ( ${freqthresh} * ${recsamplesize} )
 group by don_branch_id, don_branch_name, rec_branch_name, rec_branch_id, event_type 
) as weg
order by nb_lineages desc, avg_support desc;
""" > ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} 
wc -l ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} 
${ptgscripts}/plot_spetree_event_density.r ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} ${speciestreeBS}.lsd_internalPopulations.nwk
done &


###########################################
## 07. compare gene evolution scenarios
###########################################

export comparerecs=${entdb}/07.compare_scenarios
mkdir -p ${comparerecs}
export compoutdir=${comparerecs}/${parsedreccol}
mkdir -p ${compoutdir}/


## analyse of co-evolution!!
export minevfreqmatch=0.5
export minjoinevfreqmatch=1.0
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi

## look for correlated gene lineage histories through identification of matching speciation and transfer events

### OPTION: exclude oldest species tree branches to avoid unspecific matches (and speed-up search):
if [ ! -z maxreftreeheight ]
  # e.g.: maxreftreeheight=0.25
  exclbrlist=${coretree}/branches_older_than_${maxreftreeheight}
  python ${ptgscripts}/list_branches.py --intree ${speciestreeBS}.lsd_internalPopulations.nwk --root_age 1.0 --older_than ${maxreftreeheight} --out ${exclbrlist}
  exclbr=$(cat ${exclbrlist} | tr '\n' ',' | sed -e "s/,$/\n/g" | sed -e "s/,/', '/g")

  # create smaller table with only desired event
  sqlite3 ${sqldb} << EOF
  ALTER TABLE gene_lineage_events RENAME TO gene_lineage_events_full;
  CREATE TABLE gene_lineage_events AS  
  (SELECT gene_lineage_events_full.* 
   FROM gene_lineage_events_full 
   INNER JOIN species_tree_events USING (event_id)
   INNER JOIN species_tree on rec_branch_id=branch_id
   WHERE branch_name NOT IN ('${exclbr}') )
   AND reconciliation_id=${parsedreccolid}
  ;
  CREATE INDEX ON gene_lineage_events (replacement_label_or_cds_code);
  CREATE INDEX ON gene_lineage_events USING HASH (replacement_label_or_cds_code);
  CREATE INDEX ON gene_lineage_events (event_id);
  CREATE INDEX ON gene_lineage_events USING HASH (event_id);
  CREATE INDEX ON gene_lineage_events (freq);
  ALTER TABLE gene_lineage_events ADD PRIMARY KEY (replacement_label_or_cds_code, event_id);
  ANALYZE;
  .quit
EOF

fi
### end OPTION

# collect data
# BEWARE: GENERATES AN AWFUL LOT OF DATA, PREPARE DISK SPACE ACCORDINGLY
# indication: with defaults settings evtypeparse='ST'; minevfreqmatch=0.5; minjoinevfreqmatch=1.0; maxreftreeheight=0.25
# on a 880 Enterobacteriaceae dataset, results in ~300 GB output (made to be split into ~1GB files)
python $ptgscripts/compare_collapsedALE_scenarios.py --events_from_postgresql_db ${sqldbname} \
 --event_type ${evtypeparse} --min_freq ${minevfreqmatch} --min_joint_freq ${minjointevfreqmatch} --threads 8 \
 --dir_table_out ${compoutdir} &> $entlogs/compare_collapsedALE_scenarios.${parsedreccol}.log &

# should add constraint like reconciliation_id=${parsedreccolid} to ensure events are not matched across collections

###############################################
## 08. orthologous and clade-specific gene sets
###############################################

export orthogenes=${rapdb}/08.orthologs
mkdir -p ${orthogenes}

cd ${ptgrepo} ; export ptgversion=$(git log | grep commit | cut -d' ' -f2) ; cd -

# classify genes into orthologous groups for each gene of the reconciled gene tree sample
# do not report detailed results but, using graph analysis, combine the sample-wide classification into one classification for the gene family
# run in parallel

## for the moment only coded for dated ALE model (ALEml reconciliations)
## and with no colpased gene tree clades (need to discard events below replacement clade subtree roots [as in parse_collapsedALE_scenarios.py]
## and to transpose orthologus group classification of collapsed clade to member genes)

if [-z ${getOrpthologuesOptions} ] ; then
  getOrpthologuesOptions=" --ale.model='dated' --methods='mixed' --max.frac.extra.spe=0.5 --majrule.combine=0.5 --colour.combined.tree"
fi
if [-z ${orthoColId} ] ; then
  orthocolid=1
fi
# generate Ortholog Collection
orthocol=ortholog_collection_${orthoColId}
mkdir -p ${orthogenes}/${orthocol}
${ptgscripts}/get_orthologues_from_ALE_recs.py -i ${outrecdir} -o ${orthogenes}/${orthocol} ${getOrpthologuesOptions} -T 8 &> $raplogs/get_orthologues_from_ALE_recs_${orthocol}.log

# import ortholog classification into database
sqlite3 ${sqldb} """INSERT INTO ortholog_collections (ortholog_col_id, ortholog_col_name, reconciliation_id, software, version, algorithm, ortholog_col_date, notes) VALUES 
(${orthocolid}, '${orthocol}', ${parsedreccolid}, 'pantagruel/scripts/get_orthologues_from_ALE_recs.py', '${ptgversion:0:7}', 'getOrthologues(method=''mixed'')', '2018-07-20', 
'source from https://github.com/flass/pantagruel/commits/${ptgversion}, call: ''scripts/get_orthologues_from_ALE_recs.py ${getOrpthologuesOptions}''')
;
"""
python ${ptgscripts}/pantagruel_sqlitedb_load_orthologous_groups.py ${sqldb} ${orthogenes}/ortholog_collection_1 "mixed" "majrule_combined_0.500000" ${orthocolid}


# generate abs/pres matrix
orthocol=ortholog_collection_${orthocolid}
echo $orthocol
orthomatrad=${orthogenes}/${orthocol}/mixed_majrule_combined_0.5.orthologs
python ${ptgscripts}/get_ortholog_presenceabsence_matrix_from_sqlitedb.py ${sqldb} ${orthomatrad} ${coregenome}/${focus}/${focus}_genome_codes ${orthocolid}
${ptgscripts}/get_clade_specific_genes.r ${orthomatrad}_genome_counts.no-singletons.mat ${sqldb} ${orthocolid} ${coregenome}/${focus}/${focus} ${orthomatrad}

# list clade-specific orthologs
export orthomatrad=${orthogenes}/${orthocol}/mixed_majrule_combined_0.5.orthologs

