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

## prepare protein families for alignment
export protfamseqs=${mmseqsclout}_clusters_fasta
export protali=${ptgdb}/02.gene_alignments
nrprotali=$protali/nr_protfam_clustalo_alignments
mkdir -p ${nrprotali}/ ${ptglogs}/clustalo/
# generate lists of many jobs to be executed by a single process (clustalo jobs are too short to be dispatched via a cluster array job)
# clustalomega is already parallelized (for part of the computation, i.e. distance matrix calculation) 
# so less threads than avalilable cores (e.g. use 1 process / 2 cores) are specified
# so each process can run at least at 100% and mobilize 2 cores when parrallel,
# or more when cores unused by other concurent pocesses (during their sequential phase)
python ${ptgscripts}/schedule_ali_task.py $protfamseqs.tab $protfamseqs $protali/$(basename ${protali})_tasklist $((`nproc` / 2))

## align non-redundant protein families
for tasklist in `ls $protali/${protali##*.}_tasklist.*` ; do
  bntl=`basename $tasklist`
  ${ptgscripts}/run_clustalo_sequential.sh $tasklist $nrprotali 2 &> $ptglogs/clustalo/$bntl.clustalo.log &
done

## check consistency of non-redundant protein sets
mkdir -p $ptgtmp
protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'nr_protein_id' | cut -d':' -f1)
if [ -z $protidfield ] ; then 
 protidfield=$(head -n 1 ${genomeinfo}/assembly_info/allproteins_info.tab |  tr '\t' '\n' | grep -n 'protein_id' | cut -d':' -f1)
fi
cut -f $protidfield ${genomeinfo}/assembly_info/allproteins_info.tab | grep -v "^$\|protein_id" | sort -u > ${genomeinfo}/assembly_info/allproteins_in_gff
grep '>' ${allfaarad}.nr.faa | cut -d' ' -f1 | cut -d'>' -f2 | sort -u > ${allfaarad}.nr_protlist
# compare original dataset of nr protein (as described in the input GFF files) to the aligned nr proteome
diff ${genomeinfo}/assembly_info/allproteins_in_gff ${allfaarad}.nr_protlist > $ptgtmp/diff_prot_info_fasta
if [ -s $ptgtmp/diff_prot_info_fasta ] ; then 
  >&2 echo "ERROR $(dateprompt): inconsistent propagation of the protein dataset:"
  >&2 echo "present in aligned fasta proteome / absent in info table generated from input GFF:"
  >&2 grep '>' $ptgtmp/diff_prot_info_fasta | cut -d' ' -f2
  >&2 echo "present in info table generated from input GFF / absent in aligned fasta proteome:"
  >&2 grep '<' $ptgtmp/diff_prot_info_fasta | cut -d' ' -f2
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
  for fam in `cat $ptgtmp/pal2nal_missed_fams` ; do
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
