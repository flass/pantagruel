#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=8
export defmem=124
export defwth=24
export defhpctype='PBS'
export defchunksize=100
export withpython='true'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

# safer specifying the reconciliation collaction to parse; e.g.: 'ale_collapsed_dated_1'
if [ ! -z "$2" ] ; then
  reccol="$2"
else
  # if not inferred from the record of the last reconciliation computation
  reccol=$(cut -f4 ${alerec}/reccol)
fi

######################################################
## 07.2 Parse gene tree / Species tree reconciliations
######################################################
 
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families

### parse the inferred scenarios
# parameters to be set
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi
# derived parameters
export parsedreccol=${reccol}_parsed_${parsedreccolid}
export parsedrecs=${alerec}/parsed_recs/${parsedreccol}

outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}

mkdir -p ${parsedrecs}
tasklist=$outrecdir/${reccol}_${rectype}_uml_rec_list
${ptgscripts}/lsfullpath.py "${outrecdir}/${reccol}_${rectype}/*ml_rec" > $tasklist

parsecollogd=${ptgdb}/logs/parsecol
parsecollogs=${parsecollogd}/parse_collapsedALE_scenarios.log

ALEsourcenote=$(cut -f3 ${alerec}/reccol)


# submitted parallel job
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
[ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}

jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
rm -f ${tasklist}_${dtag}_taskchunks
for jobrange in ${jobranges[@]} ; do
  replrun="${dtag}_${jobrange}"
  tail -n +$(echo $jobrange | cut -d'-' -f1) ${tasklist} | head -n ${chunksize} > ${tasklist}_${replrun}
  echo ${tasklist}_${replrun} >> ${tasklist}_${dtag}_taskchunks
done
Nchunk=`wc -l ${tasklist}_${dtag}_taskchunks | cut -f1 -d' '`
if [ ${Nchunk} -gt 2 ] ; then
  # run as an array job
  arrayspec=" -J 1-${Nchunk}"
else
  arrayspec=""
fi

case "${hpctype}" in 
  'PBS')
    [ ${Nchunk} -gt 2 ] && arrayspec=" -J 1-${Nchunk}"
    subcmd="qsub${arrayspec} -N parseColALEscenarios -l select=1:ncpus=${ncpus}:mem=${mem}gb${parallelflags},walltime=${wth}:00:00 -o ${repllogd} -j oe -V ${ptgscripts}/parseColALEscenario_array_PBS.qsub"
    ;;
  'LSF')
    if [ ${Nchunk} -gt 2 ] ; then
      arrayspec="[1-${Nchunk}]"
	  [ ! -z ${maxatonce} ] && arrayspec="${arrayspec}%${maxatonce}"
      arrayjobtag='.%I'
    fi
    if [ ${wth} -le 12 ] ; then
      bqueue='normal'
    elif [ ${wth} -le 48 ] ; then
      bqueue='long'
    else
     bqueue='basement'
    fi
    memmb=$((${mem} * 1024)) 
    nflog="${repllogd}/parseColALEscenarios.%J${arrayjobtag}.o"
    [ -z "${parallelflags}" ] && parallelflags="span[hosts=1]"
    subcmd="bsub -J parseColALEscenarios${arrayspec} -R \"select[mem>${memmb}] rusage[mem=${memmb}] ${parallelflags}\" \
            -n${ncpus} -M${memmb} -q ${bqueue} \
            -o ${nflog} -e ${nflog} -env 'all' ${ptgscripts}/parseColALEscenarios_array_LSF.bsub"
    ;;
  *)
    echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
    exit 1
	;;
esac
echo "$subcmd"
eval "$subcmd"

#### NOTE
## here for simplicity only the log variables refering to parsed reconciliations (parsedreccol, parsedreccolid, parsedreccoldate) are recorded in the database
## but not the log variables refering to the actual reconciliations (reccol, reccolid, reccoldate)
####
export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}" > ${alerec}/parsedreccol

## store reconciliation parameters and load parsed reconciliation data into database
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh "${database}" "${sqldb}" "${parsedrecs}" "${ALEversion}" "${ALEalgo}" "${ALEsourcenote}" "${parsedreccol}" "${parsedreccolid}" "${parsedreccoldate}"
EOF
