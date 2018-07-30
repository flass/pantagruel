#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

export raproot=$1
envsourcescript=${raproot}/environ_pantagruel.sh
source $envsourcescript

collapseCladeOptions=$2

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
