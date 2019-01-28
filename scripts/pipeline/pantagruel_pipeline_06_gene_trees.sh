#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}

# to set misc variables ; not safe though
if [ -e ${ptgtmp}/nondefvardecl.sh ] ; then
  source ${ptgtmp}/nondefvardecl.sh
fi

#############################################################
## 06. Gene trees (full [ML] and collapsed [bayesian sample])
#############################################################

mkdir -p ${mlgenetrees}
mkdir -p ${ptglogs}/raxml/gene_trees

basequery="select gene_family_id, size from gene_family_sizes where gene_family_id is not null and gene_family_id!='$cdsorfanclust'"
 
 
### prepare HPC version
if [[ "$hpcremoteptgroot" != 'none' ]] ; then

  export hpcremotehost=$(echo "$hpcremoteptgroot" | cut -d':' -f1)
  export hpcremotefolder=$(echo "$hpcremoteptgroot" | cut -d':' -f2)

  python ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${protali} --base.query="${basequery}" \
   --famsets.min.sizes="4,500,2000,10000" --famsets.max.sizes="499,1999,9999,"

  # sync input and ouput folders with remote HPC host
  ${ptgscripts}/sync_ptgdb_with_remote_host.sh ${ptgdbname} ${ptgroot} \
   ${hpcremoteptgroot} ${cdsalifastacodedir} ${colalinexuscodedir}/${collapsecond} ${mlgenetrees} ${bayesgenetrees} ${protali}/cdsfams_* ${sqldb}

  echo "please connect to remote host $hpcremotehost and execute the following scripts in order "
  echo "(waiting for completion of all array jobs submitted by one script before executing the next):"
  echo "- pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh"
  echo "- pantagruel_pipeline_06.2_HPC_collapse_gene_trees.sh"
  echo "- pantagruel_pipeline_06.3_HPC_bayesian_gene_trees.sh"
  echo "- pantagruel_pipeline_06.4_HPC_replace_spe_by_pops.sh"
  echo "then copy back ouput files and updated database file by syncing the root folder from remote host to this host"


  exit 0
  
fi

### local version
    
python ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${protali} --base.query="${basequery}" \
--famsets.min.sizes="4"

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A1: no collapsing, just convert the alignments from fasta to nexus to directly compute bayesian trees
  export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code  #
  for aln in `ls ${alifastacodedir}` ; do
    convalign -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  export nexusaln4chains=${colalinexuscodedir}
  export mboutputdir=${bayesgenetrees}
  
  #### end OPTION A1
else
  #### OPTION B1: compute ML gene trees and collapse rake clades

  if [[ "${chaintype}" != 'collapsed' ]] ; then
   echo "Error: incorrect value for variable chaintype: '${chaintype}'"
   exit 1
  fi
  ##########################
  ## 06.1 Full ML gene trees
  ##########################
  
  ## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
  allcdsfam2phylo=${protali}/cdsfams_minsize4
  tasklist=${cdsfam2phylo}_aln_list
  rm -f $tasklist
  for fam in `cut -f1 $cdsfam2phylo` ; do
   ls ${cdsalifastacodedir}/${fam}.codes.aln >> $tasklist
  done
  ${ptgscripts}/raxml_sequential.sh ${tasklist} ${mlgenetrees} 'GTRCAT' 'bipartitions rootedTree identical_sequences' 'x' $(nproc) 'true'
  
  ############################
  ## 06.2 Gene tree collapsing
  ############################

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
  
  python ${ptgscripts}/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist} --diraln=${cdsalifastacodedir} --fmt_aln_in='fasta' \
   --threads=${ncpus} --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output --dir_identseq=${mlgenetrees}/identical_sequences \
   ${collapsecriteriondef}

  export collapsecoldate=$(date +%Y-%m-%d)
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
  export mboutputdir=${bayesgenetrees}/${collapsecond}
  
fi
#### end OPTION B1

############################
## 06.3 Bayesian gene trees
############################

## run mrbayes on collapsed alignments
export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
mkdir -p ${mboutputdir}
nchains=4
nruns=2
ncpus=$(( $nchains * $nruns ))
tasklist=${nexusaln4chains}_ali_list
rm -f $tasklist
${ptgscripts}/lsfullpath.py ${nexusaln4chains} > $tasklist

${ptgscripts}/mrbayes_sequential.sh ${tasklist} ${mboutputdir} 'Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'

################################################################################
## 06.4 Convert format of Bayesian gene trees and replace species by populations
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  export coltreechains=${genetrees}/${chaintype}_tree_chains
  echo "Error: not implemented yet:"
  echo "must generalize the script replace_species_by_pop_in_gene_trees.py so to only convert the format of gene tree chains, not replacing anything in them"
  #~ python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -o ${coltreechains} --threads=${ncpus} --reuse=0 --verbose=0 --logfile=${repllogs}_${replrun}.log &
  exit 1
  
  #### end OPTION A2: 
else
  #### OPTION B2: collapsed rake clades in gene trees need to be replaced by mock population leaves
  #### will edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
  if [-z ${replacecolid} ] ; then
   replacecolid=1
  fi
  export coltreechains=${genetrees}/${chaintype}_tree_chains
  export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'
  mkdir -p ${coltreechains}/${collapsecond}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  mkdir -p ${ptgdb}/logs/replspebypop
  tasklist=${coltreechains}_${collapsecond}_nexus_list
  ls $bayesgenetrees/${collapsecond}/*run1.t > $tasklist
  repllogd=${ptgdb}/logs/replspebypop
  repllogs=$repllogd/replace_species_by_pop_in_gene_trees
  replrun=$(date +'%d%m%Y')

  # local parallel run
  python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.newick -o ${coltreechains}/${collapsecond} \
   --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk --population_node_distance=${speciestree%.*}_interNodeDistPopulations \
   --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${colmethod} --threads=${ncpus} --reuse=0 --verbose=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log

  export replacecoldate=$(date +%Y-%m-%d)

  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${sqldb} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} ${collapsecriteriondef} ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}

fi
#### end OPTION B2
