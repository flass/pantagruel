#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [[[[ncpus] mem(gb)] maxwalltime(h)] paralleflags(\"str:str:...\")]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  export ncpus="$2"
else
  export ncpus=1
fi
if [ ! -z "$3" ] ; then
  mem="$3"
else
  mem=96gb
fi
if [ ! -z "$4" ] ; then
  wth="$4"
else
  wth=24
fi
if [ ! -z "$5" ] ; then
  parallelflags=":$5"
  # for instance, parallelflags="mpiprocs=1:ompthreads=8"
  # will force the use of OpenMP multi-threading instead of default MPI parallelism
  # note this is not standard and might not be the default on your HPC system;
  # also the flags may be different depending on the HPC system config.
  # you should get in touch with your system admins to know the right flags
else
  parallelflags=""
fi

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

mkdir -p ${alerec}

### perform reconciliations with ecceTERA

# parameters to be set: defaults:
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
    ls ${inputtrees}
  done > ${tasklist}
fi
teralogs=${ptgdb}/logs/ecceTERA
mkdir -p ${teralogs}/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  spetree=${speciestree}.lsd.nwk ${recsamplesize} ${ALEalgo}
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree=Stree.nwk
fi


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

if [[ -z "${terabin}" ]] ; then
  terabin="$(dirname $(readlink -f $(command -v ecceTERA)))"
fi

Njob=`wc -l $tasklist | cut -f1 -d' '`
qsubvars="tasklist='${tasklist}', resultdir='${outrecdir}', spetree='${spetree}', nrecs='${recsamplesize}', alealgo='${ALEalgo}', alebin='${alebin}', terabin='${terabin}', watchmem='${watchmem}'"
echo "qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=${ncpus}:mem=${mem}${parallelflags},walltime=${wth}:00:00 -o ${teralogs}/${reccol} -j oe -v \"$qsubvars\" ${ptgscripts}/ecceTERA_array_PBS.qsub"
qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=${ncpus}:mem=${mem}${parallelflags},walltime=${wth}:00:00 -o ${teralogs}/${reccol} -j oe -v "$qsubvars" ${ptgscripts}/ecceTERA_array_PBS.qsub

export reccoldate=$(date +%Y-%m-%d)

terasourcenote=""
pathterabin=$(readlink -f ${terabin})
terasourcenote="using ecceTERA software compiled from source; $(ecceTERA | grep version); binaries found at ${pathterabin}"
echo -e "${reccolid}\t${reccoldate}\t${terasourcenote}" > ${genetrees}/reccol

