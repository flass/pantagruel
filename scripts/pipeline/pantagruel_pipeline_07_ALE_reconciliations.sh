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
checkfoldersafe ${alerec}

if [ -z "${ptgthreads}" ] ; then
  export ptgthreads=$(nproc)
fi

###############################################
## 07. Gene tree / Species tree reconciliations
###############################################

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

### perform reconciliations with ALE

# parameters to be set: defaults:
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${ALEalgo} == 'ALEml_undated' ] ; then
  export rectype='undat'
elif  [ ${ALEalgo} == 'ALEml' ] ; then
  export rectype='dated'
else
  echo "Error: '${ALEalgo}' is not a known algorithm/executable file name for ALE software; exit now"
  exit 1
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_ALE_recs

export reccoldate=$(date +%Y-%m-%d)
if [[ -z "${alebin}" ]] ; then
  alebin=$(command -v ${ALEalgo})
  if [ -z "${alebin}" ] ; then
    echo "'${ALEalgo}' command is not available from the PATH"
    dockerale="$(docker image ls | grep alesuite | awk '{print $1,$3}')"
	if [ ! -z "${dockerale}" ] ; then
	  ALEsourcenote="using ALE Docker image ${dockerale}"
	else
	  echo "alesuite docker image not detected on this host"
	  echo "could not determine what ALE executables are being used"
	  ALEsourcenote='ALE source not determined'
	fi
  fi
fi

ALEheader="$(${ALEalgo} -h | head -n 1)"
if [ -z "${ALEheader}" ] ; then
  # will be interpreted as ALE v0.4 or less
  ALEheader="using ALE software"
fi
ALEvtag=$(echo ${ALEheader} | awk '{ print $NF }')

if [[ ! -z "${alebin}" ]] ; then
  pathalebin=$(readlink -f "${alebin}")
  alerepo=${pathalebin%%ALE/*}ALE/
  if [ -d ${alerepo} ] ; then
	alesrcvers=$(cd ${alerepo} && git log | head -n 1 | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
	alesrcorig=$(cd ${alerepo} && git remote -v | grep fetch | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  fi
  if [ ! -z "${alesrcvers}" ] ; then
    ALEsourcenote="using ALE software (version ${ALEvtag}) compiled from source; code origin: ${alesrcorig}; code version ${alesrcvers}"
  else
    ALEsourcenote="${ALEheader#*${ALEalgo} } binaries found at '${pathalebin}'"
  fi
fi

tasklist=${alerec}/${collapsecond}_${replmethod}_Gtrees_list
if [ -z ${genefamlist} ] ; then
  ${ptgscripts}/lsfullpath.py "${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk" > ${tasklist}
else
  rm -f ${tasklist}
  for fam in $(cut -f1 ${genefamlist}) ; do
    ls ${coltreechains}/${collapsecond}/${replmethod}/${fam}*-Gtrees.nwk 2> /dev/null
  done > ${tasklist} 
fi
alelogs=${ptgdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir

cd ${ptgtmp} 

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  spetree=${speciestree}.lsd.nwk
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree='Stree.nwk'
fi

if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resumetasklist
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfgs in $(cat ${tasklist}) ; do
    bng=$(basename ${nfgs})
	if [[ "${ALEvtag}" == 'v0.5' ||  "${ALEvtag}" > 'v0.5' ]] ; then
      [ ${spetree} == 'Stree.nwk' ] && aleoutSpref=${bng/Gtrees/Stree} || aleoutSpref=$(basename ${spetree})
      bnalerec=${aleoutSpref}_${bng}.ale.${tag}ml_rec
    else
	  # assume ALE v0.4 or less
	  bnalerec=${bng}.ale.${tag}ml_rec
	fi
	if [[ ! -e ${recs}/${collapsecond}/${replmethod}/${reccol}/${bnalerec} ]] ; then
     echo ${nfgs}
   fi
  done > ${tasklist}_resumetasklist
  tasklist=${tasklist}_resumetasklist
fi

## perform reconciliations sequentially (one gene family after another)
if [ -s ${tasklist} ] ; then
  export worklocal='false'
  ${ptgscripts}/ale_sequential.sh ${tasklist} ${outrecdir} ${spetree} ${recsamplesize} ${ALEalgo}
else
  echo "no gene tree left to reconcile, skip reconciliation computation"
fi

echo -e "${reccolid}\t${reccoldate}\t${ALEsourcenote}\t${reccol}" > ${alerec}/reccol
echo -e "\n# Reconciliation collection details:"
cat ${alerec}/reccol

######################################################
## 07.2 Parse gene tree / Species tree reconciliations
######################################################
#### NOTE
## here for simplicity only the log variables refering to parsed reconciliations (parsedreccol, parsedreccolid, parsedreccoldate) are recorded in the database
## but not the log variables refering to the actual reconciliations (reccol, reccolid, reccoldate)
####

### parse the inferred scenarios
# parameters to be set
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi
# derived parameters
export parsedreccol=${reccol}_parsed_${parsedreccolid}
export parsedrecs=${alerec}/parsed_recs/${parsedreccol}

mkdir -p ${parsedrecs}

reclist=${outrecdir}_rec_list
${ptgscripts}/lsfullpath.py "${outrecdir}/*ml_rec" > ${reclist}

if [ "$chaintype" == 'fullgenetree' ] ; then
  pops=""
else
  pops=" --populations ${speciestree/.full/}_populations"
fi
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families
stepa="parsing ALE scenarios"
echo "${stepa}"
python2.7 ${ptgscripts}/parse_collapsedALE_scenarios.py --rec_sample_list ${reclist} \
 ${pops} --reftree ${speciestree}.lsd.nwk --ALE_algo ${rectype} \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads ${ptgthreads}  &> ${ptglogs}/parse_collapsedALE_scenarios.log
checkexec "Could not complete ${stepa}" "Successfully ${stepa/ing/ed}"

export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}\t${parsedreccol}" > ${alerec}/parsedreccol
echo -e "\n# Parsed reconciliation collection details:"
cat ${alerec}/parsedreccol

if [ "${resumetask}" == 'true' ] ; then
  stepb="clean the database from previous inserts and indexes"
  echo "Resume mode: first ${stepb}"
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_clean_reconciliations.sh "${database}" "${sqldb}" "${parsedreccolid}"
  checkexec "Could not complete ${stepb/clean/cleaning}" "Successfully ${stepb/clean/cleaned}"
fi
stepc="Storing reconciliation parameters and load parsed reconciliation data into database"
echo ${stepc}
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh "${database}" "${sqldb}" "${parsedrecs}" "${ALEversion}" "${ALEalgo}" "${ALEsourcenote}" "${parsedreccol}" "${parsedreccolid}" "${parsedreccoldate}"
checkexec "Could not complete ${stepc}" "Successfully ${stepc/ing/ed}"

# rapid survey of event density over the reference tree
for freqthresh in 0.1 0.25 0.5 ; do
sqlite3 -cmd ".mode tabs" ${sqldb} """
select don_branch_id, don_branch_name, rec_branch_id, rec_branch_name, event_type, nb_lineages, cum_freq, cum_freq/nb_lineages as avg_support from (
 select don_branch_id, don_stree.branch_name as don_branch_name, rec_branch_id, rec_stree.branch_name as rec_branch_name, event_type, count(*) as nb_lineages, cast(sum(freq) as real)/${recsamplesize} as cum_freq
  from gene_lineage_events 
  inner join species_tree_events using (event_id) 
  inner join species_tree as rec_stree on rec_branch_id=rec_stree.branch_id
  left join species_tree as don_stree on don_branch_id=don_stree.branch_id
 where freq >= ( ${freqthresh} * ${recsamplesize} )
 group by don_branch_id, don_branch_name, rec_branch_name, rec_branch_id, event_type 
) as weg
order by nb_lineages desc, avg_support desc;
""" > ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} 
nevt=$(wc -l ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh})
echo ${nevt}
if [ $(echo ${nevt} | cut -f1 -d' ') -gt 0 ] ; then
  ${ptgscripts}/plot_spetree_event_density.r ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} ${speciestree/.full/}_collapsedPopulations.nwk
fi
done

