#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019


if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}

envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}



#############################################################
## 06.2 Full ML gene trees on HPC
#############################################################


#### OPTION: edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
  # no collapsing, just convert the alignments from fasta to nexus
  for aln in `ls ${alifastacodedir}` ; do
    convalign -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  export nexusaln4chains=${colalinexuscodedir}
  export mboutputdir=${bayesgenetrees}
  
else
 if [[ "${chaintype}" == 'collapsed' ]] ; then
  if [ -z ${collapsecolid} ] ; then
    collapsecolid=1
  fi
  if [[ "$collapseCladeParams" != 'default' ]] ; then
    eval "$collapseCladeParams"
    # e.g.:  eval 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'
  fi
  
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
  qsub -N mark_unresolved_clades -l select=1:ncpus=${ncpus}:mem=16gb,walltime=4:00:00 -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe -V -S /usr/bin/bash << EOF
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
  
 else
  echo "Error: incorrect value for variable chaintype: '${chaintype}'"
  exit 1
 fi
fi
#### end OPTION
