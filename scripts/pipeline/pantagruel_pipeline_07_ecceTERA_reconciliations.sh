#!/usr/bin/env bash

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
## 07. Gene tree / Species tree reconciliations    with ecceTERA
###############################################

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

### perform reconciliations with ecceTERA

# parameters to be set: defaults:
#~ export ALEversion='v0.4'
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${ecceTERAalgo} == 'amalgamate' ] ; then
  # work on gene tree samples
  export rectype='ecceTERA_amalgamatedGchain'
  inputtrees="${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk"
  inputtreetag="Gtrees"
else
  # work on single ML gene tree with branch supports
  inputtrees="${bayesgenetrees}/${collapsecond}/*.con.tre"
  inputtreetag="contre"
  if [ ${ecceTERAalgo:0:8} == 'collapse' ] ; then
    # ecceTERA uses branch supports to 'collapse' the unresolved clades (collapsing appraoch similar in principle but different from Pantagruel's')
    export rectype="ecceTERA_collapsedGconstree_${ecceTERAalgo##*_}"
  else
    export rectype='ecceTERA'
  fi
fi
export reccol="${chaintype}_${rectype}_${reccolid}"

tasklist=${alerec}/${collapsecond}_${replmethod}_${inputtreetag}_list
if [ -z ${genefamlist} ] ; then
  ${ptgscripts}/lsfullpath.py "${inputtrees}" > ${tasklist}
else
  for fam in $(cut -f1 ${genefamlist}) ; do
    ls ${coltreechains}/${collapsecond}/${replmethod}/${fam}*-Gtrees.nwk 2> /dev/null
  done > ${tasklist}
fi
alelogs=${ptgdb}/logs/ecceTERA
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir

if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resumetasklist
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfgs in $(cat ${tasklist}) ; do
    bng=$(basename ${nfgs})
    bnalerec=${bng}*.reconciliationsFile_canonical_symmetric.txt
    if [ ! -e ${recs}/${collapsecond}/${replmethod}/${reccol}/${bnalerec} ] ; then
     echo ${nfgs}
    fi
  done > ${tasklist}_resumetasklist
  tasklist=${tasklist}_resumetasklist
fi

cd ${ptgtmp}/

## perform receonciliations sequentially (one gene family after another)
export worklocal='false'
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  ${ptgscripts}/ecceTERA_sequential.sh ${tasklist} ${outrecdir} ${speciestree}.lsd.nwk ${ecceTERAalgo}
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  ${ptgscripts}/ecceTERA_sequential.sh ${tasklist} ${outrecdir} Stree.nwk ${ecceTERAalgo}
fi
export reccoldate=$(date +%Y-%m-%d)

if [[ -z "${terabin}" ]] ; then
  terabin="$(command -v ecceTERA)"
fi
terasourcenote=""
pathterabin=$(readlink -f ${terabin})
terasourcenote="using ecceTERA software compiled from source; $(ecceTERA | grep version); binaries found at ${pathalebin}"
echo -e "${reccolid}\t${reccoldate}\t${terasourcenote}\t${reccol}" > ${alerec}/reccol
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
${ptgscripts}/lsfullpath.py "${outrecdir}/*reconciliationsFile_canonical_symmetric.txt" > $reclist


if [ "$chaintype" == 'fullgenetree' ] ; then
  pops=""
else
  pops=" --populations ${speciestree/.full/}_populations"
fi


## below is not ready yet
echo "parsing script not coded yet; will stop here."
exit 2

## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families
python2.7 ${ptgscripts}/parse_collapsedTERA_scenarios.py --rec_sample_list ${reclist} \
 ${pops} --reftree ${speciestree}.lsd.nwk \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads ${ptgthreads}  &> ${ptglogs}/parse_collapsedTERA_scenarios.log

checkexec "Could not complete parsing ecceTERA scenarios" "Successfully parsed ecceTERA scenarios"

export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}\t${parsedreccol}" > ${alerec}/parsedreccol
echo -e "\n# Parsed reconciliation collection details:"
cat ${alerec}/parsedreccol

if [ "${resumetask}" == 'true' ] ; then
  echo "Resume mode: first clean the database from previous inserts and indexes"
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_clean_reconciliations.sh "${database}" "${sqldb}" "${parsedreccolid}"
fi
echo "Store reconciliation parameters and load parsed reconciliation data into database"
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh "${database}" "${sqldb}" "${parsedrecs}" "${ALEversion}" "${ALEalgo}" "${ALEsourcenote}" "${parsedreccol}" "${parsedreccolid}" "${parsedreccoldate}"
