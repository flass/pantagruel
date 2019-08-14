#!/bin/bash

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
checkfoldersafe ${protali}


####################################
## 02. Homologous Sequence Alignemnt
####################################

## prepare protein families for alignment
mkdir -p ${nrprotali}/ ${ptglogs}/clustalo/
tasklist=$protali/$(basename ${protfamseqs})_tasklist
python ${ptgscripts}/schedule_ali_task.py ${protfamseqs}.tab ${protfamseqs} ${tasklist} 1 "${famprefix}P000000"

## align non-redundant protein families
run_clustalo_sequential () {
  task=$1
  bn=`basename ${task}`
  fout=${bn/fasta/aln}
  echo "task: $task" &> ${ptglogs}/clustalo/${bn}.clustalo.log
  date +"%d/%m/%Y %H:%M:%S" &> ${ptglogs}/clustalo/${bn}.clustalo.log
  clustalo --threads=1 -i ${task} -o ${nrprotali}/${fout} &> ${ptglogs}/clustalo/${bn}.clustalo.log
  date +"%d/%m/%Y %H:%M:%S" &> ${ptglogs}/clustalo/${bn}.clustalo.log
}
export -f run_clustalo_sequential
parallel --joblog ${ptglogs}/run_clustalo_sequential.log run_clustalo_sequential :::: ${tasklist}

# check that alignments are not empty
${ptgscripts}/lsfullpath.py "${nrprotali}/*" > ${nrprotali}_list
rm -f ${nrprotali}_list_empty
for ali in `cat ${nrprotali}_list` ; do
  if [ ! -s $ali ] ; then
    echo "'$ali' file is empty"
    ls $ali >> ${nrprotali}_list_empty
  fi
done
if [ -e ${nrprotali}_list_empty ] ; then
  echo "Error: some alignment failed:"
  cat ${nrprotali}_list_empty
  exit 1
fi

## reconstruct full (redundant) protein alignments
# make list of CDS sets
for cg in `cat ${genomeinfo}/assemblies_list` ; do ls $cg/*_cds_from_genomic.fna.gz >> ${genomeinfo}/all_cds_fasta_list ; done

# generate (full protein alignment, unaligned CDS fasta) file pairs and reverse-translate alignments to get CDS (gene family) alignments
mkdir -p ${ptglogs}/extract_full_prot_and_cds_family_alignments/
python ${ptgscripts}/extract_full_prot_and_cds_family_alignments.py --nrprot_fam_alns ${nrprotali} --singletons ${protfamseqs}/${protorfanclust}.fasta \
 --prot_info ${genomeinfo}/assembly_info/allproteins_info.tab --repli_info ${genomeinfo}/assembly_info/allreplicons_info.tab --assemblies ${assemblies} \
 --dirout ${protali} --famprefix ${famprefix} --logs ${ptglogs}/extract_full_prot_and_cds_family_alignments --identical_prots ${allfaarad}.identicals.list
checkexec "$(promptdate)-- Critical error during the production of full CDS alignments from the nr protein alignments" "$(promptdate)-- complete generation of full CDS alignments without critical errors"

## check consistency of full reverse translated alignment set
ok=1
for fam in `ls $protali/full_cdsfam_alignments/ | cut -d'.' -f1` ; do 
nseqin1=`grep -c '>' $protali/full_cdsfam_fasta/$fam.fasta`
nseqin2=`grep -c '>' $protali/full_protfam_alignments/$fam.aln`
nseqout=`grep -c '>' $protali/full_cdsfam_alignments/$fam.aln`
if [[ "$nseqin1" != "$nseqout" || "$nseqin2" != "$nseqout" ]] ; then 
  echo -e "$fam\tfull_cdsfam_fasta: ${nseqin1}; full_protfam_alignments: ${nseqin2}; full_cdsfam_alignments: $nseqout"
  ok=0
fi
done > ${protali}/pal2nal_missed_fams
if [ $ok -lt 1 ] ; then
  >&2 promptdate 
  >&2 echo "WARNING: failure of pal2nal.pl reverse translation step for families: $(cut -f1 ${protali}/pal2nal_missed_fams | xargs)"
  >&2 echo "  (See list in ${protali}/pal2nal_missed_fams)"
  >&2 echo "  Will use tranposeAlignmentProt2CDS.py instead, a less safe, but more permissive method for generating CDS alignment"
  # some protein alignments do not match the CDS sequences
  # transpose the CDS into the positions of the aligned protein; assumes no indels, only mismatches and possibly shortenned sequences
  rm -f ${ptglogs}/tranposeAlignmentProt2CDS.log && touch ${ptglogs}/tranposeAlignmentProt2CDS.log
  for fam in `cut -f1 ${protali}/pal2nal_missed_fams` ; do
    ${ptgscripts}/tranposeAlignmentProt2CDS.py $protali/full_cdsfam_fasta/$fam.fasta $protali/full_protfam_alignments/$fam.aln $protali/full_cdsfam_alignments/$fam.aln > ${ptglogs}/tranposeAlignmentProt2CDS.log
  done
  checkexec "$(promptdate)-- failed to generate the reverse translated aligments missed by pal2nal" "$(promptdate)-- Complete generating the reverse translated aligments missed by pal2nal"
  >&2 promptdate 
  >&2 echo "See '${ptglogs}/tranposeAlignmentProt2CDS.log' for the list of alignments that were reverse-translated using the coarse algorithm implemented in tranposeAlignmentProt2CDS.py"
fi

# join non-ORFan and ORFan family count matrices
rm -f ${protali}/all_families_genome_counts.mat*
cat ${protali}/full_families_genome_counts-noORFans.mat > ${protali}/all_families_genome_counts.mat
tail -n +2 ${protali}/${famprefix}C000000_genome_counts-ORFans.mat >> ${protali}/all_families_genome_counts.mat
gzip ${protali}/all_families_genome_counts.mat
checkexec "$(promptdate)-- failed to generate the gene family count matrices" "$(promptdate)-- Complete generating the gene family count matrices"

