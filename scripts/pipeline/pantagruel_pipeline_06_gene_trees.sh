#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [gene_fam_list]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  export genefamlist="$2"
fi

# function for checking the success of every step
checkexec (){
  if [ $? != 0 ]; then
    echo "ERROR: $1" 1>&2
    exit 1
  else
    if [ ! -z "$2" ] ; then
      echo "$2"
    fi
  fi
}


#############################################################
## 06. Gene trees (full [ML] and collapsed [bayesian sample])
#############################################################

mkdir -p ${colalinexuscodedir}/${collapsecond}/ ${mlgenetrees}/ ${bayesgenetrees}/
mkdir -p ${ptglogs}/raxml/gene_trees

basequery="select gene_family_id, size from gene_family_sizes where gene_family_id is not null and gene_family_id!='${cdsorfanclust}'"
 
 
### prepare HPC version
if [[ ! -z "$hpcremoteptgroot" && "$hpcremoteptgroot" != 'none' ]] ; then

  export hpcremotehost=$(echo "$hpcremoteptgroot" | cut -d':' -f1)
  export hpcremotefolder=$(echo "$hpcremoteptgroot" | cut -d':' -f2)
  
  # try and size the job regarding to the gene tree sizes: do separate lists by gene family size
  python ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${genetrees} --base.query="${basequery}" \
   --famsets.min.sizes="4,500,2000,10000" --famsets.max.sizes="499,1999,9999,"
  allfamlists=$(ls ${genetrees}/cdsfams_*)
  checkexec "could not generate gene family lists ; exit now" "succesfully generated gene family lists : $allfamlists"
  
  if [ -z "$genefamlist" ] ; then
    famlists="${allfamlists}"
  else
    echo "will restrict the gene families to be processed to those listed in '$genefamlist':"
    cut -f1 ${genefamlist} | head
    ngf=$(wc -l ${genefamlist} | cut -d' ' -f1)
    if [ $ngf -gt 6 ] ; then
      echo "... ($ngf families)"
    fi
    bngenefamlist=$(basename ${genefamlist})
    rm -f ${genetrees}/${bngenefamlist}_*
    for fam in $(cut -f1 ${genefamlist}) ; do
      for cdsfamlist in ${allfamlists} ; do
        grep ${fam} ${cdsfamlist} >>  ${genetrees}/${bngenefamlist}_${cdsfamlist}
      done
    done
    famlists=$(ls ${genetrees}/${bngenefamlist}_cdsfams_*)
    checkexec "could not generate restricted gene family lists ; exit now" "succesfully generated restricted gene family lists : $famlists"
  fi
  
  # sync input and ouput folders with remote HPC host
  ${ptgscripts}/sync_ptgdb_with_remote_host.sh ${ptgdbname} ${ptgroot} \
   ${hpcremoteptgroot} ${cdsalifastacodedir} ${colalinexuscodedir}/${collapsecond} ${mlgenetrees} ${bayesgenetrees} ${famlists} ${sqldb}

  echo "please connect to remote host $hpcremotehost and execute the following scripts in order "
  echo "(waiting for completion of all array jobs submitted by one script before executing the next):"
  echo "- pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh"
  echo "- pantagruel_pipeline_06.2_HPC_collapse_gene_trees.sh"
  echo "- pantagruel_pipeline_06.3_HPC_bayesian_gene_trees.sh"
  echo "- pantagruel_pipeline_06.4_HPC_replace_spe_by_pops.sh"
  echo "- pantagruel_pipeline_06.5_HPC_populate_db_collapsed_clades.sh"
  echo "then copy back ouput files and updated database file by syncing the root folder from remote host to this host"

  exit 0
  
fi

### local version
    
allfamlist=${genetrees}/cdsfams_minsize4
python ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${genetrees} --base.query="${basequery}" \
 --famsets.min.sizes="4"
checkexec "could not generate gene family list ; exit now" "succesfully generated gene family list : $allfamlist"

if [ -z "$genefamlist" ] ; then
  famlist=${allfamlist}
else
  echo "will restrict the gene families to be processed to those listed in '$genefamlist':"
  cut -f1 ${genefamlist} | head
  ngf=$(wc -l ${genefamlist} | cut -d' ' -f1)
  if [ $ngf -gt 6 ] ; then
    echo "... ($ngf families)"
  fi
  bngenefamlist=$(basename ${genefamlist})
  famlist=${genetrees}/${bngenefamlist}_cdsfams_minsize4
  rm -f ${famlist}
  for fam in $(cut -f1 ${genefamlist}) ; do
    grep ${fam} ${allfamlist} >> ${famlist}
  done
  checkexec "could not generate restricted gene family list ; exit now" "succesfully generated restricted gene family lists : $famlist"
fi

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A1: no collapsing, just convert the alignments from fasta to nexus to directly compute bayesian trees
  for aln in `ls ${alifastacodedir}` ; do
    convalign -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  checkexec "could not convert alignemts from Fasta to Nexus format ; exit now" "succesfully converted alignemts from Fasta to Nexus format"
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  checkexec "could not move Nexus alignemts to destination folder ${colalinexuscodedir}/ ; exit now" "succesfully moved Nexus alignemts to destination folder ${colalinexuscodedir}/"
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
  tasklist=${famlist}_aln_list
  rm -f ${tasklist}
  for fam in $(cut -f1 ${famlist}) ; do
   ls ${cdsalifastacodedir}/${fam}.codes.aln >> ${tasklist}
  done
  ${ptgscripts}/raxml_sequential.sh ${tasklist} ${mlgenetrees} 'GTRCATX' 'bipartitions rootedTree identical_sequences' 'x' $(nproc) 'true'
  checkexec "RAxML tree estimation was interupted ; exit now" "RAxML tree estimation complete"
  
  ############################
  ## 06.2 Gene tree collapsing
  ############################

  if [ -z ${collapsecolid} ] ; then
    collapsecolid=1
  fi
  
  ## detect clades to be collapsed in gene trees
  export collapsecriteriondef="--clade_stem_conds=[('$criterion','>=',$cladesupp)] --within_clade_conds=[('$withinfun','$criterion','<=',$subcladesupp,-1),('max','$criterion','<',$cladesupp,-1)]"
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  mlgenetreelist=${mlgenetrees%*s}_list
  ${ptgscripts}/lsfullpath.py "${mlgenetrees}/${mainresulttag}/*" | sort > ${mlgenetreelist}
  
  python ${ptgscripts}/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist} --diraln=${cdsalifastacodedir} --fmt_aln_in='fasta' \
   --threads=$(nproc) --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output --dir_identseq=${mlgenetrees}/identical_sequences \
   ${collapsecriteriondef} 
  checkexec "ML tree collapsing was interupted ; exit now" "ML tree collapsing complete"

  export collapsecoldate=$(date +%Y-%m-%d)
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
  export mboutputdir=${bayesgenetrees}/${collapsecond}
  echo -e "${collapsecolid}\t${collapsecoldate}\t${collapsecriteriondef}" > ${genetrees}/replacecol
fi
#### end OPTION B1

############################
## 06.3 Bayesian gene trees
############################

## run mrbayes on collapsed alignments
#~ export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
mkdir -p ${mboutputdir}
nchains=4
nruns=2
ncpus=$(( $nchains * $nruns ))
tasklist=${nexusaln4chains}_ali_list
rm -f $tasklist
${ptgscripts}/lsfullpath.py "${nexusaln4chains}/*" > $tasklist

${ptgscripts}/mrbayes_sequential.sh ${tasklist} ${mboutputdir} 'Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'
checkexec "MrBayes tree estimation was interupted ; exit now" "MrBayes tree estimation complete"

################################################################################
## 06.4 Convert format of Bayesian gene trees and replace species by populations
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  echo "Error: not implemented yet:"
  echo "must generalize the script replace_species_by_pop_in_gene_trees.py so to only convert the format of gene tree chains, not replacing anything in them"
  #~ python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -o ${coltreechains} --threads=${ncpus} --reuse=0 --verbose=0 --logfile=${repllogs}_${replrun}.log &
  exit 1
  
  #### end OPTION A2: 
else
  #### OPTION B2: collapsed rake clades in gene trees need to be replaced by mock population leaves
  #### will edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
  if [ -z ${replacecolid} ] ; then
   replacecolid=1
  fi
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
   --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${colmethod} --threads=$(nproc) --reuse=0 --verbose=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log
  checkexec "replacement of collapsed clades was interupted ; exit now" "replacement of collapsed clades complete"

  export replacecoldate=$(date +%Y-%m-%d)
  echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol
  
  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${sqldb} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} "${collapsecriteriondef}" ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}
  checkexec "populating the SQL database with collapsed/replaced gene tree clades was interupted ; exit now" "populating the SQL database with collapsed/replaced gene tree clades complete"

fi
#### end OPTION B2
