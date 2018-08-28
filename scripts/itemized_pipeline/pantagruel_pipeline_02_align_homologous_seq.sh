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

####################################
## 02. Homologous Sequence Alignemnt
####################################

export protfamseqs=${mmseqsclout}_clusters_fasta
export protali=${ptgdb}/02.gene_alignments
## prepare protein families for alignment
nrprotali=$protali/nr_protfam_clustalo_alignments
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
parrallel --joblog ${ptglogs}/run_clustalo_sequential.log run_clustalo_sequential :::: ${tasklist}

## reconstruct full (redundant) protein alignments
# make list of CDS sets
for cg in `cat ${genomeinfo}/assemblies_list` ; do ls $cg/*_cds_from_genomic.fna.gz >> ${genomeinfo}/all_cds_fasta_list ; done

# generate (full protein alignment, unaligned CDS fasta) file pairs and reverse-translate alignments to get CDS (gene family) alignments
mkdir -p ${ptglogs}/extract_full_prot_and_cds_family_alignments/
python ${ptgscripts}/extract_full_prot_and_cds_family_alignments.py --nrprot_fam_alns ${nrprotali} --singletons ${protfamseqs}/${protorfanclust}.fasta \
 --prot_info ${genomeinfo}/assembly_info/allproteins_info.tab --repli_info ${genomeinfo}/assembly_info/allreplicons_info.tab --assemblies ${assemblies} \
 --dirout ${protali} --famprefix ${famprefix} --logs ${ptglogs}/extract_full_prot_and_cds_family_alignments --identical_prots ${allfaarad}.identicals.list

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
done > $ptgtmp/pal2nal_missed_fams
if [ $ok -lt 1 ] ; then
  >&2 echo "WARNING $(dateprompt): failure of pal2nal.pl reverse translation step for families: $(cat $ptgtmp/pal2nal_missed_fams | xargs)"
  >&2 echo "will use tranposeAlignmentProt2CDS.py instead, a less safe, but more permissive method for generating CDS alignment"
  # some protein alignments do not match the CDS sequences
  # transpose the CDS into the positions of the aligned protein; assumes no indels, only mismatches and possibly shortenned sequences
  for fam in `cat $ptgtmp/pal2nal_missed_fams` ; do
    ${ptgscripts}/tranposeAlignmentProt2CDS.py $protali/full_cdsfam_fasta/$fam.fasta $protali/full_protfam_alignments/$fam.aln $protali/full_cdsfam_alignments/$fam.aln
  done
fi

# join non-ORFan and ORFan family count matrices
rm -f ${protali}/all_families_genome_counts.mat*
cat ${protali}/full_families_genome_counts-noORFans.mat > ${protali}/all_families_genome_counts.mat
tail -n +2 ${protali}/${famprefix}C000000_genome_counts-ORFans.mat >> ${protali}/all_families_genome_counts.mat
gzip ${protali}/all_families_genome_counts.mat
