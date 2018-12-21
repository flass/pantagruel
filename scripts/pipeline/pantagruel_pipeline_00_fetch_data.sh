#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$2" ] ; echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"
export ptgroot="$2"
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}

#################################
## 00. Data download and grooming
#################################

if [ ! -e ${ncbitax}/names.dmp ] ; then
echo "did not find the relevant taxonomy flat files in '${ncbitax}/'; dowlod the from NCBI Taxonomy FTP"
  ### Download data
  ## using FTP
  # downloading relevant tables from NCBI Taxonomy database 
  user='anonymous'
  pswd=$myemail
  host='ftp.ncbi.nlm.nih.gov'
  openparam=" -u $user,$pswd $host"
  source="/pub/taxonomy"
  files="taxcat.tar.gz* taxcat_readme.txt taxdump.tar.gz* taxdump_readme.txt"
  mkdir -p $ncbitax
  lftp -c "open $openparam; cd $source ; mget -O $ncbitax/ $files ; quit"
  # reduce taxonomic id2name reference file complexity
  for tgz in `ls $ncbitax/*.tar.gz` ; do md5sum -c $tgz.md5 && tar -xzf $tgz ; done
fi

## using NCBI web interface:
# Search Assembly database using a query defined as convenient:	
# e.g.: 'txid543[Organism:exp] AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter]) AND ("complete genome"[filter])'
# save assemblies (Download Assemblies > Source Database: RefSeq; File type: All file type (including assembly-structure directory)) as $indata/genome_assemblies_DATASETTAG.tar
# extract assembly data
mkdir -p  $indata/assemblies/
cd $ncbiass/
for tarf in `ls genome_assemblies_*.tar` ; do 
  tartag=${tarf#genome_assemblies_*}
  datasettag=${tartag%.*}
  tar -tf genome_assemblies_${datasettag}.tar | cut -d'/' -f2 | sort -u | grep -v "README\|report" > $ncbiass/genome_assemblies_${datasettag}_list
  tar -xf genome_assemblies_${datasettag}.tar 
  mv report.txt report_${datasettag}.txt
  assd=(`ls ncbi-genomes-* -d`)
  for dass in `ls ${assd[0]}/ | grep -v README` ; do
   if [ ! -z $dass ] ; then
    rm -rf $ncbiass/$dass
    mv ${assd[0]}/$dass $ncbiass/
   fi
  done
  rm -r ${assd[0]}/ 
done
grep "# Organism name" $ncbiass/*/*_assembly_stats.txt > $ncbiass/all_assemblies_organism_names

# store in centralised folder for NCBI assemblies and just record links
mkdir -p $indata/assemblies/
for ass in `cat $ncbiass/genome_assemblies_*_list` ; do ln -s $ncbiass/${ass}* $indata/assemblies/  ; done


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
prokka-genbank_to_fasta_db ${indata}/assemblies/*/*gbff > ${ptgtmp}/${refgenus}.faa >& ${ptgtmp}/prokka-genbank_to_fasta_db.log
# cluster similar sequences
cdhit -i ${ptgtmp}/${refgenus}.faa -o ${ptgtmp}/${refgenus}_representative.faa -T 0 -M 0 -G 1 -s 0.8 -c 0.9 &> cdhit.log
# 805428 representative proteins
rm -fv ${ptgtmp}/${refgenus}.faa ${ptgtmp}/${refgenus}_representative.faa.clstr
# replace database name for genus detection by Prokka 
prokkablastdb=$(dirname $(dirname $(ls -l `which prokka` | awk '{ print $NF }')))/db/genus/
cp -p ${ptgtmp}/${refgenus}_representative.faa ${prokkablastdb}/${refgenus}
cd ${prokkablastdb}/
makeblastdb -dbtype prot -in ${refgenus}
cd -

# run Prokka 
export annot=${customassemb}/prokka_annotation
export straininfo=${customassemb}/strain_infos.txt
mkdir -p $annot
for allcontigs in `ls ${customassemb}/contigs/` ; do
gproject=${allcontigs%%_*}
date
echo "### assembly: $gproject; contig files from: ${customassemb}/contigs/${allcontigs}"
echo "running Prokka..."
$ptgscripts/run_prokka.sh ${gproject} ${customassemb}/contigs/${allcontigs} ${straininfo} ${annot}/${gproject} ${seqcentre}
echo "fix prokka output to integrate region information into GFF files"
annotoutgff=($(ls ${annot}/${gproject}/*.gff))
python $ptgscripts/add_region_feature2prokkaGFF.py ${annotoutgff[0]} ${annotoutgff[0]%*.gff}.ptg.gff ${straininfo} ${customassemb}/contigs/${allcontigs} ${assembler}
echo "fix prokka output to integrate taxid information into GBK files"
annotoutgbk=($(ls ${annot}/${gproject}/*.gbk))
python $ptgscripts/add_taxid_feature2prokkaGBK.py ${annotoutgbk[0]} ${annotoutgbk[0]%*.gbk}.ptg.gbk ${straininfo}
echo "done."
date
done &> ${ptgtmp}/${ptgdbname}_customassembly_annot_prokka.log &

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
mkdir -p ${manuin}/
touch ${manuin}/manual_metadata_dictionary.tab ${manuin}/manual_curated_metadata_dictionary.tab ${manuin}/manual_dbxrefs.tab
rm -rf ${genomeinfo}/assemblies*
ls ${assemblies}/* -d > ${genomeinfo}/assemblies_list
mkdir -p ${genomeinfo}/assembly_metadata
## extract assembly/sample metadata from flat files
python ${ptgscripts}/extract_metadata_from_gbff.py --assembly_folder_list=${genomeinfo}/assemblies_list --add_raw_metadata=${manuin}/manual_metadata_dictionary.tab \
--add_curated_metadata=${manuin}/manual_curated_metadata_dictionary.tab --add_dbxref=${manuin}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism" --output=${genomeinfo}/assembly_metadata

export ngenomes=$((`wc -l ${genomeinfo}/assembly_metadata/metadata.tab | cut -d' ' -f1` - 1))

mv ${ptgdb}/environ_pantagruel_${ptgdbname}.sh ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0 && \
 sed -e "s#'REPLACEngenomes'#$ngenomes#" ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0 > ${ptgdb}/environ_pantagruel_${ptgdbname}.sh && \
 rm ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0

echo "work with a database of $ngenomes genomes (excluding lab strains)"
