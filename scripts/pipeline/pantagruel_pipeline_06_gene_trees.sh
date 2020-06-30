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

if [ -z "${ptgthreads}" ] ; then
  export ptgthreads=$(nproc)
fi


if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  echo "WARNING: these environment variables are set: chaintype='${chaintype}' ; recmethod='${recmethod}'."
  echo " -> running this task (06/gene_trees) of Pantagruel pipeline when using GeneRax reconciliation method is DEPRECATED."
  echo ""
  echo "If you really want to run task 06 under these settings, please set the env. variable 'raxmlfullgenetrees' as 'true' with:"
  echo "   \`export raxmlfullgenetrees=true\`"
  echo "(but note that these gene trees will not be used by the next task)"
  echo ""
  echo "Finally, you may consider using the -c|--collapse option of the Pantagruel pipeline,"
  echo "which allows to reconstruct the history of large pangenomes with high sequence diversity and complex histories,"
  echo "in which case task 06 is required to pre-compute gene trees so to collapse their rake clades."
  echo ""
  echo -e "\n# # #\n"
  if [ "${raxmlfullgenetrees}" != 'true' ] ; then
    echo ""
    echo "This task will exit with code 0 now, moving on to the next task (07/reconciliations)."
    exit 0
  else
    echo "raxmlfullgenetrees='true' ; Will execute task 06 regardless."
  fi
fi

###############################################
## 06. Gene trees (full and collapsed ML trees)
###############################################

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
  echo "- pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh [OPTIONS] ${hpcremoteptgroot}/$(basename ${envsourcescript})"
  echo "- pantagruel_pipeline_06.2_HPC_collapse_gene_trees.sh [OPTIONS] ${hpcremoteptgroot}/$(basename ${envsourcescript})"
  echo "- pantagruel_pipeline_06.3_HPC_bayesian_gene_trees.sh [OPTIONS] ${hpcremoteptgroot}/$(basename ${envsourcescript})"
  echo "- pantagruel_pipeline_06.4_HPC_replace_spe_by_pops.sh [OPTIONS] ${hpcremoteptgroot}/$(basename ${envsourcescript})"
  echo "- pantagruel_pipeline_06.5_HPC_populate_db_collapsed_clades.sh [OPTIONS] ${hpcremoteptgroot}/$(basename ${envsourcescript})"
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

##########################
## 06.1 Full ML gene trees
##########################

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
tasklist=${famlist}_aln_list
rm -f ${tasklist}*
for fam in $(cut -f1 ${famlist}) ; do
  aln=${cdsalifastacodedir}/${fam}.codes.aln
  if [[ "${resumetask}" == "true" ]] ; then
    mainres=${mlgenetrees}/${mainresulttag}/RAxML_${mainresulttag}.${fam}.codes
    skiptag=${mlgenetrees}/bulk/${fam}.codes.smallreducedali
    if [[ ! -s ${mainres} && ! -s ${skiptag} ]] ; then
      ls ${aln} >> ${tasklist}_resume
    fi
  fi
  ls ${aln} >> ${tasklist}
done
step1="step 1: RAxML estimation of ML gene trees"
step1s="step 1: gene tree estimation"
echo ${step1}
if [[ "${resumetask}" == "true" ]] ; then
  if [[ -s ${tasklist}_resume ]] ; then
    echo "Resume mode: $(wc -l ${tasklist}_resume | cut -d' ' -f1) ML trees left to infer"
  else
    echo -e "Resume mode: all ML tree built; skip ${step1s}\n# # #"
  fi
  tasklist=${tasklist}_resume
else
  echo "step 1: build ML gene trees"
fi
if [ -s "${tasklist}" ] ; then
  ${ptgscripts}/raxml_sequential.sh "${tasklist}" "${mlgenetrees}" 'GTRCATX' 'bipartitions rootedTree identical_sequences' 'x' "${ptgthreads}" 'true'
  checkexec "${step1s} was interupted. You can check the cause in logs stored in '${ptglogs}/raxml/gene_trees/' ; exit now" "${step1s} complete\n# # #"
fi
  
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A1: no collapsing, just convert the alignments from fasta to nexus to directly compute bayesian trees

  echo -e "step 2: env var \$chaintype is set to '${chaintype}', hence will not collapse the gene trees; nothing to do here.\n# # #"
  
  #### end OPTION A1
else
  #### OPTION B1: collapse rake clades

  if [[ "${chaintype}" != 'collapsed' ]] ; then
   echo "Error: incorrect value for variable chaintype: '${chaintype}'"
   exit 1
  fi

  ############################
  ## 06.2 Gene tree collapsing
  ############################
  step2="step 2: collapse rake clades in ML gene trees"
  step2s="step 2: ML gene tree collapsing"
  echo ${step2}
  if [ -z ${collapsecolid} ] ; then
    if [ -e ${genetrees}/collapsecol ] ; then
      collapsecolid=$(cut -f1 ${genetrees}/collapsecol)
      if [[ "${resumetask}" != "true" ]] ; then
        collapsecolid=$(( ${collapsecolid} + 1 ))
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
    rm -f ${mlgenetreelist}_resume
    for gt in $(cat ${mlgenetreelist}) ; do
      bngt=$(basename $gt)
      fam=$(echo ${bngt} | cut -d'.' -f2 )
      colaln=${colalinexuscodedir}/${collapsecond}/collapsed_alns/${fam}-collapsed.nex
      colmlgt=${colalinexuscodedir}/${collapsecond}/collapsed_ML_genetrees/${fam}-collapsed.nwk
      if [[ ! -s ${colaln} || ! -s ${colmlgt} ]] ; then
        echo ${gt} >> ${mlgenetreelist}_resume
      fi
    done
    if [[ -s ${mlgenetreelist}_resume ]] ; then
      echo "Resume mode: $(wc -l ${mlgenetreelist}_resume | cut -d' ' -f1) ML trees left to collapse"
    else
      echo -e "Resume mode: all ML tree collapsed; skip ${step2s}\n# # #"
    fi
    mlgenetreelist=${mlgenetreelist}_resume
  fi

  if [ -s "${mlgenetreelist}" ] ; then
    python2.7 ${ptgscripts}/mark_unresolved_clades.py --in_gene_tree_list=${mlgenetreelist} --diraln=${cdsalifastacodedir} --fmt_aln_in='fasta' \
     --threads=${ptgthreads} --dirout=${colalinexuscodedir}/${collapsecond} --no_constrained_clade_subalns_output --dir_identseq=${mlgenetrees}/identical_sequences \
     ${collapsecriteriondef} 
    checkexec "${step2s} was interupted ; exit now" "${step2s} complete\n# # #"
  fi
  export collapsecoldate=$(date +%Y-%m-%d)
#  export mboutputdir=${bayesgenetrees}/${collapsecond}
  echo -e "${collapsecolid}\t${collapsecoldate}\t${collapsecriteriondef}" > ${genetrees}/collapsecol
fi
#### end OPTION B1

################################################################################
## 06.3 Replace species by populations
################################################################################

# backward compatability fix:
if [[ -d ${genetrees}/${chaintype}_tree_chains && ! -d ${coltreechains} ]] ; then
  cd ${genetrees}/
  ln -s ${chaintype}_tree_chains $(basename ${coltreechains})
  cd - > /dev/null
else
  mkdir -p ${coltreechains}/
fi

repltasklist=${genetrees}/$(basename ${colmlgenetrees})_list
${ptgscripts}/lsfullpath.py ${colmlgenetrees}/ > ${repltasklist}
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  ## no need to replace anything in the tree; just link the original ML gene trees to the folder of replaced gene tree chains
  echo "step 3: env var \$chaintype is set to '${chaintype}', hence no need for gene tree format conversion and replacement of collapsed clades"
  echo -e "will directly use ML trees from RAxML; link files from '${colmlgenetrees}/' into '${coltreechains}/${collapsecond}/${replmethod}/'\n# # #"
  rm -rf ${coltreechains}/${collapsecond}/${replmethod}/
  mkdir -p ${coltreechains}/${collapsecond}/${replmethod}/
  for nfgt in $(cat ${repltasklist}) ; do
    bnfgt=$(basename ${nfgt})
    radnfgt=${bnfgt#*.}  # trim RAxML_tag. prefix
    ln -s ${nfgt} ${coltreechains}/${collapsecond}/${replmethod}/${radnfgt}-Gtree.nwk
  done
  #### end OPTION A2
else
  #### OPTION B2: collapsed rake clades in gene trees need to be replaced by mock population leaves
  step3="step 3: replace collapsed clades in ML gene trees by population representatives"
  step3s="step 3: replacement of collapsed clades"
  echo ${step3}
  repllogd=${ptgdb}/logs/replspebypop
  mkdir -p ${repllogd}/
  repllogs=${repllogd}/replace_species_by_pop_in_gene_trees
  replrun=$(date +'%d%m%Y')  
  if [ "${resumetask}" == 'true' ] ; then
    # resume mode (useful fter a stop in batch computing, or to collect those jobs that crashed and may need to be re-ran with more mem/time allowance)
	# evaluate what gene tree parsing/replacement jobs remain to be done
    rm -f ${repltasklist}_resume
    for nfcolgt in $(cat ${repltasklist}) ; do
      bncolgt=$(basename ${nfcolgt})
      bnGtre=${bncolgt/.nwk/-replaced-Gtree.nwk}
      bnGaln=${bncolgt/.nwk/-replaced.aln}
      bnStre=${bncolgt/.nwk/-Stree.nwk}
      if [[ ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnGtre} || ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnGaln} || ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnStre} ]] ; then
        echo ${nfcolgt}
      fi
    done > ${repltasklist}_resume
    if [ -s ${repltasklist}_resume ] ; then
      echo "Resume step 3: $(wc -l ${repltasklist}_resume | cut -d' ' -f1) collapsed ML trees remain to be processed for replacement of collapsed clades"
    else
      echo "Resume step 3: all collapsed ML trees processed; skip ${step3s}\n# # #"
    fi
    repltasklist=${repltasklist}_resume
  fi
  #### will edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history  
  if [ -z ${replacecolid} ] ; then
    if [ -e ${genetrees}/replacecol ] ; then
      replacecolid=$(cut -f1 ${genetrees}/replacecol)
      if [[ "${resumetask}" != "true" ]] ; then
        replacecolid=$(( ${replacecolid} + 1 ))
      fi
    else
      replacecolid=1
    fi
  fi
  if [ -s ${repltasklist} ] ; then
    mkdir -p ${coltreechains}/${collapsecond}
    ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
    # local parallel run
    python2.7 ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${repltasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.nwk \
	 -o ${coltreechains}/${collapsecond} --flatRCs --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk \
	 --population_node_distance=${speciestree%.*}_interNodeDistPopulations --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${replmethod} \
	 --threads=${ptgthreads} --reuse=0 --verbose=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log
    checkexec "${step3s} was interupted ; exit now" "${step3s} complete\n# # #"
  fi 
  ## make a summary matrix of collapsed clade (CC) occurence in genome populations
  ## this really reflects the output of step 2, but relies on files made during step 4
  matCC=${coltreechains}/${collapsecond}/${replmethod}_PopFreqsCC.mat
  echo -e -n "\t" > ${matCC}
  grep -v '#' ${coretreerad}_populations | cut -f1 | tr '\n' '\t' >> ${matCC}
  echo "" >> ${matCC}
  cat ${coltreechains}/${collapsecond}/${replmethod}_phyloprofiles/* >> ${matCC}
  
  export replacecoldate=$(date +%Y-%m-%d)
  echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol

  ##############################################################
  ## 06.4 Populate database with all collapsed gene tree results
  ##############################################################
  if [ "${resumetask}" == 'true' ] ; then
    # first clean the database
    ${ptgscripts}/pantagruel_sqlitedb_phylogeny_clean_collapsed_clades.sh "${sqldb}" "${collapsecolid}" "${replacecolid}"
  fi
  step4="step 4: populating the SQL database with collapsed/replaced gene tree clades"
  step4s="step 4: population of the SQL database"
  echo ${step4}
  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh "${database}" "${sqldb}" "${colalinexuscodedir}" "${coltreechains}" "${collapsecond}" "${replmethod}" "${collapsecriteriondef}" "${collapsecolid}" "${replacecolid}" "${collapsecoldate}" "${replacecoldate}"
  checkexec "${step4s} was interupted ; exit now" "${step4s} complete\n# # #"

  #### end OPTION B2
 
fi

