#!/bin/bash

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
checkfoldersafe ${alerec}

###############################################
## 07. Gene tree / Species tree reconciliations
###############################################

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

### perform reconciliations with ALE

# parameters to be set: defaults:
#~ export ALEversion='v0.4'
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${ALEalgo} == 'ALEml_undated' ] ; then
  export rectype='undat'
else
  export rectype='dated'
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_recs

tasklist=${alerec}/${collapsecond}_${replmethod}_Gtrees_list
ls ${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk > $tasklist
alelogs=${ptgdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir

cd ${ptgtmp} 

if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resumetasklist
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfgs in $(cat $tasklist) ; do
    bng=$(basename $nfgs)
    bnalerec=${bng}.ale.${tag}ml_rec
    if [[ ! -e ${recs}/${collapsecond}/${replmethod}/${reccol}/${bnalerec} ]] ; then
     echo ${nfgs}
   fi
  done > ${tasklist}_resumetasklist
  tasklist=${tasklist}_resumetasklist
fi

## perform receonciliations sequentially (one gene family after another)
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  ${ptgscripts}/ale_sequential.sh ${tasklist} ${outrecdir} ${speciestree}.lsd.nwk ${recsamplesize} ${ALEalgo}
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  ${ptgscripts}/ale_sequential.sh ${tasklist} ${outrecdir} Stree.nwk ${recsamplesize} ${ALEalgo}
fi
export reccoldate=$(date +%Y-%m-%d)
if [ ! -z "$alebin" ] ; then
  alerepo=${alebin%%ALE/*}ALE/
  alesrcvers=$(cd ${alerepo} && git log | head -n 1 | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  alesrcorig=$(cd ${alerepo} && git remote -v | grep fetch | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  if [ ! -z "$alesrcvers" ] ; then
    ALEsourcenote="using ALE software compiled from source; code origin: ${alesrcorig}; version ${alesrcvers}"
  else
    ALEsourcenote="using ALE software binaries found at $(readlink -f ${alebin})"
  fi
else
  ALEsourcenote="using ALE Docker image $(docker image ls | grep alesuite | awk '{print $1,$3}')"
fi
echo -e "${reccolid}\t${reccoldate}\t${ALEsourcenote}" > ${genetrees}/reccol

### parse the inferred scenarios
# parameters to be set
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi
# derived parameters
export parsedreccol=${reccol}_parsed_${parsedreccolid}
export parsedrecs=${alerec}/parsed_recs/${parsedreccol}

mkdir -p ${parsedrecs}
reclist=$outrecdir/ale_collapsed_${rectype}_uml_rec_list
${ptgscripts}/lsfullpath.py "${outrecdir}/ale_collapsed_${rectype}/*ml_rec" > $reclist

if [ "$chaintype" == 'fullgenetree' ] ; then
  pops=""
else
  pops=" --populations ${speciestree/.full/}_populations"
fi
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families
python ${ptgscripts}/parse_collapsedALE_scenarios.py --rec_sample_list ${reclist} \
 ${pops} --reftree ${speciestree}.lsd.nwk \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads 8  &> $entlogs/parse_collapsedALE_scenarios.log &

export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}" > ${genetrees}/parsedreccol



######################################################
## 07.2 Parse gene tree / Species tree reconciliations
######################################################
#### NOTE
## here for simplicity only the log variables refering to parsed reconciliations (parsedreccol, parsedreccolid, parsedreccoldate) are recorded in the database
## but not the log variables refering to the actual reconciliations (reccol, reccolid, reccoldate)
####

## store reconciliation parameters and load parsed reconciliation data into database
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh ${database} ${sqldb} ${parsedrecs} ${ALEversion} ${ALEalgo} ${ALEsourcenote} ${parsedreccol} ${parsedreccolid} ${parsedreccoldate}

# rapid survey of event density over the reference tree
for freqthresh in 0.1 0.25 0.5 ; do
sqlite3 ${sqldb} """
.mode tabs 
select don_branch_id, don_branch_name, rec_branch_id, rec_branch_name, event_type, nb_lineages, cum_freq, cum_freq/nb_lineages as avg_support from (
 select don_branch_id, don_stree.branch_name as don_branch_name, rec_branch_id, rec_stree.branch_name as rec_branch_name, event_type, count(*) as nb_lineages, sum(freq)::real/${nsample} as cum_freq
  from gene_lineage_events 
  inner join species_tree_events using (event_id) 
  inner join species_tree as rec_stree on rec_branch_id=rec_stree.branch_id
  left join species_tree as don_stree on don_branch_id=don_stree.branch_id
 where freq >= ( ${freqthresh} * ${recsamplesize} )
 group by don_branch_id, don_branch_name, rec_branch_name, rec_branch_id, event_type 
) as weg
order by nb_lineages desc, avg_support desc;
""" > ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} 
wc -l ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} 
${ptgscripts}/plot_spetree_event_density.r ${parsedrecs}/summary_gene_tree_events_minfreq${freqthresh} ${speciestree/.full/}_collapsedPopulations.nwk
done &

