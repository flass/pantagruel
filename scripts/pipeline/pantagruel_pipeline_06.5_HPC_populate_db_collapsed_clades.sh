#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019


# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=1
export defmem=32
export defwth=12
export defhpctype='PBS'
export defchunksize=1
export withpython='true'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

################################################################################
## 06.5 Populate database with all collapsed gene tree results
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: pretty much nothing to do
  echo "Warning: Pantagruel task 06, step 5 ('06.5 Populate database with all collapsed gene tree results') is not required when \${chaintype}='fullgenetree' ; exit now."
  exit 0
  
  #### end OPTION A2: 
else
  source ${ptgscripts}/load_python2.7_env.sh
  # The command above will attempt to load the correct environment for Python 2.7 scripts, 
  # given what options are available on the local (HPC) system.
  # The preferred option is to use anaconda to build your own environment
  # where you make sure you are using anaconda2 with python2.7 and have all required packages installed with it.
  # You can do so using the following command:
  # conda create -n env_python2 python=2.7 python-igraph biopython bcbio-gff scipy
  
  #### OPTION B2: rake clades in gene trees were collapsed and later replaced by mock population leaves
  #### will feed data relative to these operation to the SQL databse

  ## make a summary matrix of collapsed clade (CC) occurence in genome populations [this belongs to task 06.4, but must be executed after its completion]
  matCC=${coltreechains}/${collapsecond}/${replmethod}_PopFreqsCC.mat
  echo -e -n "\t" > ${matCC}
  grep -v '#' ${coretreerad}_populations | cut -f1 | tr '\n' '\t' >> ${matCC}
  echo "" >> ${matCC}
  cat ${coltreechains}/${collapsecond}/${replmethod}_phyloprofiles/* >> ${matCC}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  export collapsecolid=$(cut -f1 ${genetrees}/collapsecol)
  export collapsecoldate=$(cut -f2 ${genetrees}/collapsecol)
  export collapsecriteriondef="$(cut -f3 ${genetrees}/collapsecol)"
  export replacecolid=$(cut -f1 ${genetrees}/replacecol)
  export replacecoldate=$(cut -f2 ${genetrees}/replacecol)
  
  echo "collapsecolid: '${collapsecolid}'; collapsecoldate: '${collapsecoldate}'; collapsecriteriondef: '${collapsecriteriondef}'; replacecolid: '${replacecolid}'; replacecoldate: '${replacecoldate}'; "
	
  ## load these information into the database
  
  case "${hpctype}" in 
   'PBS')
    subcmdpref="qsub${arrayspec} -N ptgSqlitePhyloPopulateCCs -l select=1:ncpus=${ncpus}:mem=${mem}gb${parallelflags},walltime=${wth}:00:00 -o ${repllogd} -j oe -V"
    ;;
   'LSF')
    if [ ${wth} -le 12 ] ; then
      bqueue='normal'
    elif [ ${wth} -le 48 ] ; then
      bqueue='long'
    else
     bqueue='basement'
    fi
    memmb=$((${mem} * 1024)) 
    nflog="${repllogd}/ptgSqlitePhyloPopulateCCs.%J.o"
    [ -z "${parallelflags}" ] && parallelflags="span[hosts=1]"
    subcmdpref="bsub -J replSpePopinGs${arrayspec} -R \"select[mem>${memmb}] rusage[mem=${memmb}] ${parallelflags}\" \
            -n${ncpus} -M${memmb} -q ${bqueue} -o ${nflog} -e ${nflog} -env 'all'"
    ;;
   *)
    echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
    exit 1
	;;
  esac

  subcmd = "${subcmdpref} ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh"
  echo "${subcmd}"
  eval "${subcmd}"
fi
#### end OPTION B2
