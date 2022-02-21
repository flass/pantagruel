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
checkfoldersafe ${protali}

if [ -z "${ptgthreads}" ] ; then
  ptgthreads=1
fi

####################################
## 02. Homologous Sequence Alignemnt
####################################

echo "This script relies on GNU Parallel for multi-threading"
parallel --citation

## prepare protein families for alignment
mkdir -p ${nrprotali}/ ${ptglogs}/clustalo/
tasklist=$protali/$(basename ${protfamseqs})_tasklist
python2.7 ${ptgscripts}/schedule_ali_task.py ${protfamseqs}.tab ${protfamseqs} ${tasklist} 1 "${famprefix}P000000"

## align non-redundant protein families
run_clustalo_sequential () {
  task=$1
  bn=`basename ${task}`
  fout=${bn/fasta/aln}
  echo "task: $task" &> ${ptglogs}/clustalo/${bn}.clustalo.log
  if [[ "${resumetask}" == "true" && -s ${nrprotali}/${fout} ]] ; then
	echo "output file exists; skip task"
  else
    date +"%d/%m/%Y %H:%M:%S" &> ${ptglogs}/clustalo/${bn}.clustalo.log
    clustalo --threads=1 -i ${task} -o ${nrprotali}/${fout} &> ${ptglogs}/clustalo/${bn}.clustalo.log
    date +"%d/%m/%Y %H:%M:%S" &> ${ptglogs}/clustalo/${bn}.clustalo.log
  fi
}
export -f run_clustalo_sequential
step1="step 1: alignment of nr protein families with clustalo and GNU parallel"
echo "$(promptdate)-- ${step1}"
parallel --jobs ${ptgthreads} --joblog ${ptglogs}/run_clustalo_sequential.log run_clustalo_sequential :::: ${tasklist}
checkexec "something went worng during ${step1}" "succefully completed ${step1}"

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
step20="step 2.0: generation of full CDS alignments from the nr protein alignments"
echo "$(promptdate)-- ${step2}"
python2.7 ${ptgscripts}/extract_full_prot_and_cds_family_alignments.py --nrprot_fam_alns ${nrprotali} --singletons ${protfamseqs}/${protorfanclust}.fasta \
 --prot_info ${genomeinfo}/assembly_info/allproteins_info.tab --repli_info ${genomeinfo}/assembly_info/allreplicons_info.tab --assemblies ${assemblies} \
 --dirout ${protali} --famprefix ${famprefix} --identical_prots ${allfaarad}.identicals.list \
 --threads ${ptgthreads} --logs ${ptglogs}/extract_full_prot_and_cds_family_alignments
checkexec "Critical error occurred during ${step20}" "completed ${step20} without critical errors"

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
  step21="step 2.1: generation of the reverse translated aligments missed by pal2nal"
  echo "$(promptdate)-- ${step21}"
  for fam in `cut -f1 ${protali}/pal2nal_missed_fams` ; do
    ${ptgscripts}/tranposeAlignmentProt2CDS.py $protali/full_cdsfam_fasta/$fam.fasta $protali/full_protfam_alignments/$fam.aln $protali/full_cdsfam_alignments/$fam.aln > ${ptglogs}/tranposeAlignmentProt2CDS.log
  done
  checkexec "failed during ${step21}" "Completed ${step21}"
  >&2 promptdate 
  >&2 echo "See '${ptglogs}/tranposeAlignmentProt2CDS.log' for the list of alignments that were reverse-translated using the coarse algorithm implemented in tranposeAlignmentProt2CDS.py"
fi

# join non-ORFan and ORFan family count matrices
step3="step 3: generating the gene family count matrices"
echo "$(promptdate)-- ${step4}"
rm -f ${protali}/all_families_genome_counts.mat*
cat ${protali}/full_families_genome_counts-noORFans.mat > ${protali}/all_families_genome_counts.mat
tail -n +2 ${protali}/${famprefix}C000000_genome_counts-ORFans.mat >> ${protali}/all_families_genome_counts.mat
gzip ${protali}/all_families_genome_counts.mat
checkexec "failed during ${step3}" "Completed ${step3}"

 if [[ "${compress}" == 'on' ]] ; then
   step4="compressing redundant folders"
   echo "Pantagruel compress option (-z) ON: $step4"
   cd ${protali}/
   for daln in nr_protfam_clustalo_alignments full_cdsfam_fasta ; do
     tar -czf ${daln}.tar.gz ${daln} && rm -r ${daln} || echo "Warning: could not succesfully compress '${protali}/${daln}/'; keep the full folder as is"
   done
 fi