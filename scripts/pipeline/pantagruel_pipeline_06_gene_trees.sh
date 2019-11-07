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

checkptgversion
checkfoldersafe ${genetrees}


#############################################################
## 06. Gene trees (full [ML] and collapsed [bayesian sample])
#############################################################

mkdir -p ${colalinexuscodedir}/${collapsecond}/ ${mlgenetrees}/ ${bayesgenetrees}/
mkdir -p ${ptglogs}/raxml/gene_trees

basequery="select gene_family_id, size from gene_family_sizes where gene_family_id is not null and gene_family_id!='${cdsorfanclust}'"
 
 
### prepare HPC version
if [[ ! -z "$hpcremoteptgroot" && "$hpcremoteptgroot" != 'none' ]] ; then

  if [[ "${chaintype}" == 'fullgenetree' ]] ; then
    echo "Error: HPC support not implemented yet for non-collapsed gene trees"
    echo "the computation shuld be fairly light though (only bayesian gene trees to compute) so this can be run locally"
    exit 1
  fi
  
  export hpcremotehost=$(echo "$hpcremoteptgroot" | cut -d':' -f1)
  export hpcremotefolder=$(echo "$hpcremoteptgroot" | cut -d':' -f2)
  
  # try and size the job regarding to the gene tree sizes: do separate lists by gene family size
  python2.7 ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${genetrees} --base.query="${basequery}" \
   --famsets.min.sizes="4,500,2000,10000" --famsets.max.sizes="499,1999,9999,"
  allfamlists="$(ls ${genetrees}/cdsfams_*)"
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
        grep ${fam} ${cdsfamlist} >> ${genetrees}/${bngenefamlist}_${cdsfamlist}
      done
    done
    famlists="$(ls ${genetrees}/${bngenefamlist}_cdsfams_*)"
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
python2.7 ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${genetrees} --base.query="${basequery}" \
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
  for fam in $(cut -f1 ${genefamlist}) ; do
    grep ${fam} ${allfamlist}
  done > ${famlist}
  nrestrfam=$(wc -l $famlist)
  checkexec "could not generate restricted gene family list ; exit now" "succesfully generated restricted gene family lists :\n$famlist (containing $(echo $nrestrfam | cut -d' ' -f1) gene families)"
fi

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A1: no collapsing, just convert the alignments from fasta to nexus to directly compute bayesian trees
  cdsalifastacodealnlist=${genetrees}/cdsalifastacode_aln_list
  rm -f ${cdsalifastacodealnlist}
  for fam in $(cut -f1 ${famlist}) ; do
   ls ${cdsalifastacodedir}/${fam}.codes.aln >> ${cdsalifastacodealnlist}
  done
  if [[ "${resumetask}" == "true" && -d ${colalinexuscodedir}/${collapsecond} ]] ; then
    rm -f ${cdsalifastacodealnlist}_resume && touch ${cdsalifastacodealnlist}_resume
    for aln in $(cat ${cdsalifastacodealnlist}) ; do
      bnaln=$(basename $aln)
      if [[ ! -s ${colalinexuscodedir}/${collapsecond}/${bnaln/.codes.aln/.codes.nex} ]] ; then
        echo $aln >> ${cdsalifastacodealnlist}_resume
      fi
    done
    if [[ -s ${cdsalifastacodealnlist}_resume ]] ; then
      echo "Resume task 6: $(wc -l ${cdsalifastacodealnlist}_resume | cut -d' ' -f1) alignments left to convert"
    else
      echo "Resume task 6: skip converting alignments"
    fi
    cdsalifastacodealnlist=${cdsalifastacodealnlist}_resume
  fi
  
  for aln in $(cat ${cdsalifastacodealnlist}) ; do
    bnaln=$(basename $aln)
    cp -p $aln ${colalinexuscodedir}/${collapsecond}/ && \
     convalign -i fasta -e nex -t dna nexus ${colalinexuscodedir}/${collapsecond}/$bnaln && \
      rm ${colalinexuscodedir}/${collapsecond}/$bnaln
  done
  checkexec "could not convert alignments from Fasta to Nexus format ; exit now" "succesfully converted alignemts from Fasta to Nexus format"
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}
  export mboutputdir=${bayesgenetrees}/${collapsecond}
  
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
  rm -f ${tasklist}*
  for fam in $(cut -f1 ${famlist}) ; do
    aln=${cdsalifastacodedir}/${fam}.codes.aln
    if [[ "${resumetask}" == "true" ]] ; then
      if [[ ! -s ${mlgenetrees}/${mainresulttag}/RAxML_${mainresulttag}.${fam}.codes ]] ; then
        ls ${aln} >> ${tasklist}_resume
      fi
    fi
    ls ${aln} >> ${tasklist}
  done
  if [[ "${resumetask}" == "true" ]] ; then
    if [[ -s ${tasklist}_resume ]] ; then
      echo "Resume task 6: $(wc -l ${tasklist}_resume | cut -d' ' -f1) ML trees left to infer"
    else
      echo "Resume task 6: all ML tree built; skip ML tree building"
    fi
    tasklist=${tasklist}_resume
  fi
  if [ -z "${ptgthreads}" ] ; then
    raxthreads=${ptgthreads}
  else
    raxthreads=$(nproc)
  fi
  ${ptgscripts}/raxml_sequential.sh ${tasklist} ${mlgenetrees} 'GTRCATX' 'bipartitions rootedTree identical_sequences' 'x' ${raxthreads} 'true'
  checkexec "RAxML tree estimation was interupted ; exit now" "RAxML tree estimation complete"
  
  ############################
  ## 06.2 Gene tree collapsing
  ############################
  if [ -z ${collapsecolid} ] ; then
    if [ -e ${genetrees}/replacecol ] ; then
      collapsecolid=$(cut -f1 ${genetrees}/replacecol)
      if [[ "${resumetask}" != "true" ]] ; then
        collapsecolid=$(( $collapsecolid + 1 ))
      fi
    else
      collapsecolid=1
    fi
  fi
  
  ## detect clades to be collapsed in gene trees
  export collapsecriteriondef="--clade_stem_conds=[('$criterion','>=',$cladesupp)] --within_clade_conds=[('$withinfun','$criterion','<=',$subcladesupp,-1),('max','$criterion','<',$cladesupp,-1)]"
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  mlgenetreelist=${mlgenetrees%*s}_list
  rm -f ${mlgenetreelist}
  ${ptgscripts}/lsfullpath.py "${mlgenetrees}/${mainresulttag}/*" | sort > ${mlgenetreelist}
  if [[ "${resumetask}" == "true" ]] ; then
    for gt in $(cat ${mlgenetreelist}) ; do
      bngt=$(basename $gt)
      fam=$( echo ${gt} | cut -d'.' -f2 )
      colaln=${colalinexuscodedir}/${collapsecond}/collapsed_alns/${fam}-collapsed.nex
      if [[ ! -s ${colaln} ]] ; then
        echo ${gt} >> ${mlgenetreelist}_resume
      fi
    done
    if [[ -s ${mlgenetreelist}_resume ]] ; then
      echo "Resume task 6: $(wc -l ${mlgenetreelist}_resume | cut -d' ' -f1) ML trees left to collapse"
    else
      echo "Resume task 6: all ML tree collapsed; skip ML tree collapsing"
    fi
    mlgenetreelist=${mlgenetreelist}_resume
  fi

  python2.7 ${ptgscripts}/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist} --diraln=${cdsalifastacodedir} --fmt_aln_in='fasta' \
   --threads=$(nproc) --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output --dir_identseq=${mlgenetrees}/identical_sequences \
   ${collapsecriteriondef} 
  checkexec "ML tree collapsing was interupted ; exit now" "ML tree collapsing complete"

  export collapsecoldate=$(date +%Y-%m-%d)
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
  export mboutputdir=${bayesgenetrees}/${collapsecond}
  echo -e "${collapsecolid}\t${collapsecoldate}\t${collapsecriteriondef}" > ${genetrees}/collapsecol
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
ngen=2000000
samplef=500
ncpus=$(( ${nchains} * ${nruns} ))
ntreeperchain=$(( ${ngen} / ${samplef} ))
mbtasklist=${nexusaln4chains}_ali_list
${ptgscripts}/lsfullpath.py "${nexusaln4chains}/*.nex" > ${mbtasklist}

# determine the set of numbered gene family prefixes to make separate folders
# and breakdown the load of files per folder
awk -F'/' '{print $NF}' ${mbtasklist} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}" | sort -u > ${nexusaln4chains}_ali_numprefixes
for pref in  `cat ${nexusaln4chains}_ali_numprefixes` ; do
  mkdir -p ${mboutputdir}/${pref}/
done

if [[ "${resumetask}" == "true" ]] ; then
  rm -f ${mbtasklist}_alreadydone
  rm -f ${mbtasklist}_resume
  for nfaln in $(cat ${mbtasklist}) ; do
    chaindone=''
    chainstarted=''
    nfrad1=$(basename ${nfaln})
    nfrad2=${nfrad1%.*}
    pref=$(echo ${nfrad2} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}")
    gtchain1=${mboutputdir}/${pref}/${nfrad2}.mb.run1.t
    if [[ -s ${gtchain1} ]] ; then
      if [[ $(grep -F -c "tree gen" ${gtchain1} | cut -d' ' -f1) -ge ${ntreeperchain} && ! -z "$(tail -n 1 ${gtchain1} | grep 'end')" ]] ; then
        chaindone='yes'
      else
        chainstarted='yes'
      fi
    fi
    if [ -z "${chaindone}" ] ; then
      echo ${nfaln} >> ${mbtasklist}_resume
    else
      echo ${nfaln} >> ${mbtasklist}_alreadydone
    fi
  done
fi

if [[ "${resumetask}" == "true" && -e ${mbtasklist}_alreadydone ]] ; then
  echo "$(wc -l ${mbtasklist}_alreadydone | cut -d' ' -f1) bayesian tree chains already complete; skip their computation"
  mbtasklist=${mbtasklist}_resume
  apyes='append=yes'
fi

if [ -s ${mbtasklist} ] ; then
  mbopts="Nruns=${nruns} Ngen=${ngen} Nchains=${nchains} Samplefreq=${samplef} ${apyes}"
  echo "Will now run MrBayes in parallel (i.e. sequentially for each gene alignment, with several alignments processed in parallel"
  echo "with options: ${mbopts}"
  echo ""
  ${ptgscripts}/mrbayes_sequential.sh ${mbtasklist} ${mboutputdir} "${mbopts}"
  checkexec "MrBayes tree estimation was interupted ; exit now" "MrBayes tree estimation complete"
else
  echo "no bayesian gene tree left to compute; skip to next step"
fi

################################################################################
## 06.4 Convert format of Bayesian gene trees and replace species by populations
################################################################################

mkdir -p ${ptgdb}/logs/replspebypop
repltasklist=${coltreechains}_${collapsecond}_nexus_list
${ptgscripts}/lsfullpath.py "${bayesgenetrees}/${collapsecond}/${famprefix}*/*.mb.run1.t" > ${repltasklist}
repllogd=${ptgdb}/logs/replspebypop
repllogs=$repllogd/replace_species_by_pop_in_gene_trees
replrun=$(date +'%d%m%Y')  

export dtag="$(date +'%Y%m%d-%H%M%S')"
if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resume
  # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfrun1t in $(cat $tasklist) ; do
    bnrun1t=$(basename $nfrun1t)
    bnGtre=${bnrun1t/.mb.run1.t/-Gtrees.nwk}
    if [[ "${chaintype}" == 'fullgenetree' ]] ; then
      if [ ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnGtre} ] ; then
        echo ${nfrun1t}
      fi
    else
      bnStre=${bnrun1t/.mb.run1.t/-Stree.nwk}
      if [[ ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnGtre} || ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnStre} ]] ; then
        echo ${nfrun1t}
      fi
    fi
  done > ${repltasklist}_resume
  repltasklist=${repltasklist}_resume
fi
  
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  ## no need to replace anything in the tree, just convert format from Nexus to Newick treee chains
  python2.7 ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${repltasklist} --no_replace -o ${coltreechains}/${collapsecond} --threads=${ncpus} --reuse=0 --verbose=0 --logfile=${repllogs}_${replrun}.log &
  checkexec "conversion of gene tree chains was interupted ; exit now" "conversion of gene tree chains complete"
 
  #### end OPTION A2: 
else
  #### OPTION B2: collapsed rake clades in gene trees need to be replaced by mock population leaves
  #### will edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
  if [ -z ${replacecolid} ] ; then
   replacecolid=1
  fi
  mkdir -p ${coltreechains}/${collapsecond}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone

  # local parallel run
  python2.7 ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${repltasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.nwk -o ${coltreechains}/${collapsecond} \
   --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk --population_node_distance=${speciestree%.*}_interNodeDistPopulations \
   --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${replmethod} --threads=$(nproc) --reuse=0 --verbose=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log
  checkexec "replacement of collapsed clades was interupted ; exit now" "replacement of collapsed clades complete"

  export replacecoldate=$(date +%Y-%m-%d)
  echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol

  ##############################################################
  ## 06.5 Populate database with all collapsed gene tree results
  ##############################################################
  
  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh "${database}" "${sqldb}" "${colalinexuscodedir}" "${coltreechains}" "${collapsecond}" "${replmethod}" "${collapsecriteriondef}" "${collapsecolid}" "${replacecolid}" "${collapsecoldate}" "${replacecoldate}"
  checkexec "populating the SQL database with collapsed/replaced gene tree clades was interupted ; exit now" "populating the SQL database with collapsed/replaced gene tree clades complete"

  #### end OPTION B2
fi
