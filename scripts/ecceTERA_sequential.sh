#!/bin/bash

tasklist=${1}
resultdir=${2}
spetree=${3}
teraalgo=${4}

echo "eceTERRA_sequential.sh call was: ${@}"
echo "with variables set as: tasklist=${tasklist}, resultdir=${resultdir}, spetree=${spetree}, teraalgo=${teraalgo}"

usage() {
  echo "Usage: eceTERRA_sequential.sh task_list result_dir species_tree [ecceTERA_algorithm; default='amalgamate']"
}

nbtaxafromtree () {
python << EOF
with open('${1}', 'r') as fchain:
  chainone = fchain.readline()
  print 'ntaxa:', chainone.count('(') + 2

EOF
}

## verify key variable definition
# tasklist
echo "tasklist:"
if [ -z "$tasklist" ] ; then
  echo "ERROR: need to define variable tasklist ; exit now"
  usage
  exit 2
else
  ls $tasklist
  if [ $? != 0 ] ; then
    echo "ERROR: file '$tasklist' is missing ; exit now"
  usage
    exit 2
  fi
fi

# resultdir
echo "resultdir:"
if [ -z "$resultdir" ] ; then
  echo "ERROR: need to define variable resultdir ; exit now"
  usage
  exit 2
else
  ls $resultdir -d
  if [ $? != 0 ] ; then
    echo "directory '$resultdir' is missing ; create it now"
    mkdir -p $resultdir
    if [ $? != 0 ] ; then
      echo "could not create directory '$resultdir' ; exit now"
      exit 2
    fi
  fi
fi

# spetree
echo "spetree:"
if [ -z "$spetree" ] ; then
  echo "ERROR: need to define variable spetree ; exit now"
  usage
  exit 2
fi

# nrecs
echo "nrecs:"
if [ -z "$nrecs" ] ; then
  echo -n "Default: "
  nrecs=1000
fi
echo "will sample $nrecs reconciliations"

# teraalgo
echo "teraalgo:"
if [ -z "$teraalgo" ] ; then
  echo -n "Default: "
  teraalgo='amalgamate'
fi
echo "will use $teraalgo algorithm for reconciliation estimation"

# relburninfrac
echo "relburninfrac:"
if [ -z "$relburninfrac" ] ; then
  echo -n "Default: "
  relburninfrac=0.25
fi
echo "will discard $relburninfrac fraction of the tree chain as burn-in prior to amalgamation"

# alebin (facultative location for ALE executables; default to those found in $PATH, then the Docker container)
if [ ! -z "$alebin" ] ; then
  if [ ! -z "$(ls -d "$alebin" 2> /dev/null)" ] ; then
    alebin="${alebin%*/}/"
  fi
else
  if [ -z "$(command -v ALEobserve)" ] ; then
    # no ALE program available in the command line environment (i.e. listed in $PATH directories)
    # define alebin prefix as the Docker container
    alebin="docker run -v $PWD:$PWD -w $PWD boussau/alesuite "
    # when using the Docker cntainer with local mounting, files will need to be available locally
    echo "using a Docker container with local mount, must set worklocal='true'"
    worklocal='yes'
  fi
fi
# maxpcmem
echo "# maxpcmem:"
if [ -z "${maxpcmem}" ] ; then
  maxpcmem=90
fi
echo "Will interupt reconciliation tasks if over ${maxpcmem} of server memory capacity"


# worklocal
echo "# worklocal:"
# worklocal='yes' indicates that there will be a specific machine where the 'compute work'
# is to be done 'locally' (e.g. if using a HPC system, a worker node) 
# that is different from the machine where the master/submission script is executed (e.g. login node),
# implying file traficking between them at the begin and end of the job.
if [ -z "$worklocal" ] ; then
  echo "(Use default)"
  worklocal="no"
else
  if [[ "$worklocal" == "n" || "$worklocal" == "false" ]] ; then
    worklocal="no"
  elif [[ "$worklocal" == "y" || "$worklocal" == "true" ]] ; then
    worklocal="yes"
  fi
fi
echo "will work (read/write) locally: ${worklocal}"
echo ""
echo "# # # #"
echo ""

## define ecceTERA algorithm and settings

makerunopt (){
	if [ ! -z ${4} ] ; then
		op="${4}."
	else
		op="$(basename ${2})."
	fi
	echo "species.file=${1} gene.file=${2} output.dir=${3} output.prefix=${op}"
}

commonopt="print.newick=1 print.reconciliations=1"

collapseopts=""
if [[ "${teraalgo:0:8}"=='collapse' ]] ; then
  # for single-gene tree input mode of ecceTERA:
  # 'collapse' nodes under a given threshold: the most parsimonious scenario will be considered 
  # among those inferred from all the possible binary topologies once unsupported nodes are collapsed
  # !!! this is an option of ecceTERA; it is related in principle to the collapsing performed by pantagruel (pipelin option -c) task 06, but is not equivalent
  # it is highly recomended NOT to combine those collpasing options.
  # ecceTERA built-in collpasing should be preferred,
  # unless the user desires to reduce the species tree dimension and/or to separate recent from encient gene flow,
  # in which case the pantagruel collpsing is to be preferred.
  # 
  # the collapsing threshold must be given through the $teraalgo variable: 
  #   collapse_k-xx means collpase mode is k (should always be 1 for downstream pipeline compatibility), with threshold xx. e.g. teraalgo='collapse_1-0.5'
  colopts=($(echo ${teraalgo} | cut -d'_' -f2 | tr '-' ' '))
  [ ! -z "${colopts[0]}" ] && colmode=${colopts[0]} || colmode=1
  [ ! -z "${colopts[1]}" ] && colthre=${colopts[1]} || colthre=0.5
  collapseopts="collapse.mode=${colmode} collapse.threshold=${colthre}"
fi


## loop over gene families
for nfchain in $(cat $tasklist) ; do
  echo "current task:"
  echo $nfchain
  nbtaxafromtree $nfchain
  echo ""
  echo "# # # #"
  dnchain=`dirname $nfchain`
  bnchain=`basename $nfchain`
  nfrad=${bnchain%%-*}

  ls $spetree
  if [ $? != 0 ] ; then
    echo "look for $spetree species tree file in $dnchain/ folder"
    ls ${dnchain}/${nfrad}*${spetree}*
    if [ $? != 0 ] ; then
      echo "ERROR: file '$spetree' is missing ; exit now"
      exit 2
    else
      echo "found it!" 
      lnfstree=(`ls ${dnchain}/${nfrad}*${spetree}*`)
      nfstree=${lnfstree[0]}
      echo "will use nfstree=${nfstree}"
    fi
  else
    nfstree=${spetree}
  fi
  bnstree=`basename $nfstree`

  ####
  if [[ "$worklocal"=="yes" ]] ; then
    # copy input files locally
    rsync -az ${nfchain} ${nfstree} ./
    ls -lh ${bnchain} ${bnstree}
    if [ $? != 0 ] ; then
      echo "ERROR: could not copy input file ${bnchain} or ${bnstree} locally; exit now"
      exit 2
    else
      chain="./${bnchain}"
      stree="./${bnstree}"
    fi
    # will copy output files into output dir
    savecmd="rsync -az"
  else
    chain=${nfchain}
    stree=${nfstree}
    # will rapatriate output files into output dir
    savecmd="mv -f"
  fi


  ## TERA reconciliation

  # start timing in seconds
  SECONDS=0
 
  if [[ "${teraalgo:0:10}" == 'amalgamate' ]] ; then
    # take multiple Newick gene trees as input, i.e a gene tree chain
	# this is the recomended option.
	#
	# this is a more memory-intensive option compared to single-gene tree options below,
	# but it remains relatively cheap compared to ALE probabilistic reconciliation
	
    if [[ -e ${nfchain}.ale ]] ; then
      # will use ALEobserve to produce .ale file of amalgamated gene trees from the input gene tree chain, with defined burn-in fraction
      if [[ "$worklocal" == "yes" ]] ; then
        # copy input files locally
        rsync -az ${nfchain}.ale ./
      fi
      echo "use pre-existing ALE index file:"
      ls ${nfchain}.ale
    elif [[ -e ${resultdir}/${bnchain}.ale ]] ; then
      if [[ "$worklocal" == "yes" ]] ; then
        # copy input files locally
        rsync -az ${resultdir}/${bnchain}.ale ./
	  else
	    ln -s ${resultdir}/${bnchain}.ale ${chain}.ale
      fi
      echo "use pre-existing ALE index file:"
      ls -lh ${chain}.ale
    else
      # prepare ALE index
      lenchain=`wc -l ${chain} | cut -d' ' -f1`
      burnin=`python -c "print int(${lenchain} * ${relburninfrac})"`
      echo "input tree chain is ${lenchain} long; burnin is set to ${burnin%%.*}"
      echo "# ${alebin}ALEobserve ${chain} burnin=${burnin%%.*}"
      ${alebin}ALEobserve ${chain} burnin=${burnin%%.*}
    fi
    promptdate
  
    runopt="$(makerunopt ${stree} ${chain}.ale ${resultdir})"
    teracmd="ecceTERA ${runopt} ${commonopt} ale=1 amalgamate=1"

  
  else
    # take single Newick gene tree as input, e.g. the gene tree chain consensus, or a ML tree
	if [ ! -s ${chain}.nwk ] ; then
	  # assume the input gene tree is a Nexus-formated consensus gene tree, as obtained from Mr Bayes (with 2 tree blocks)
	  # converts it from Nexus to Newick
	  python ${ptgscripts}/convert_mrbayes_constree_nex2nwk.py ${chain}
	fi
	runopt="$(makerunopt ${stree} ${chain}.nwk ${resultdir})"
	teracmd="ecceTERA ${runopt} ${commonopt} ${collapseopts}"
  
  fi
  
  echo "# ${teracmd}"
  # run it in bg with a daemon checking
  ${teracmd} &
  terapid=$!
  runmin=0
  top -b -n 1 -p ${terapid} | tail -n 2 > ${nfrad}.ecceTERA.toplog
  while [ ! -z $(ps -q ${terapid} -o comm=) ] ; do
	# fine grained record of what's happening, storing just the last value of time and mem
    TERAMEM=$(pmap ${terapid} | tail -n1 | awk '{print $NF}')
    echo "$nfrad\t$teraalgo\t$TERAMEM\tkB" > ${nfrad}.ecceTERA.memusage
    TERATIME=$SECONDS
    echo -e "$nfrad\t$teraalgo\t$TERATIME\ts" > ${nfrad}.ecceTERA.computetime
	if [ $(( $SECONDS / 60 )) -gt ${runmin} ] ; then
	  # more thorough report, logged every minute
      top -b -n 1 -p ${alepid} | tail -n 1 >> ${nfrad}.ecceTERA.toplog
	  # and sync of potential results (mostly the .ale.computetime, .ale.memusage and .ale.toplog files, as results are only written aththe end)
      ${savecmd} ./${nfrad}* ${resultdir}/
	  runmin=$(( $SECONDS / 60 ))
    fi
	pcmem=$(ps -o pid,%mem | grep ${terapid} | awk '{print $NF}')
    ## check memory use is not going off the charts
    if [ ${pcmem%.*} -ge ${maxpcmem} ] ; then
      # stop immediately
      echo "!!! Memory use is > ${pcmem%.*}% the server capacity; stop the ${nfrad} job now"
      kill -9 ${terapid}
    fi
    sleep 2s
  done
  
 
 
 
  echo ""
  echo "# # # #"

  TERATIME=$SECONDS
  if [ ! -z "${prevcomputetime}" ] ; then TERATIME=$(( ${TERATIME} + ${prevcomputetime} )) ; fi
  echo -e "${nfrad}\t${alealgo}\t${TERATIME}" > $nfrad.ecceTERA.computetime
  echo "reconciliation estimation took" $(date -u -d @${TERATIME} +"%Hh%Mm%Ss") "total time"
  if [ ! -z "${prevcomputetime}" ] ; then echo "(including ${prevcomputetime} in previous run)" ; fi

  echo "# ls ./${nfrad}*"
  ls ./${nfrad}*

  # save files
  ls ./${nfrad}* > /dev/null
  if [ ${?} == 0 ] ; then
    savecmd1="${savecmd} ./${nfrad}* ${resultdir}/"
    echo "# ${savecmd1}"
    ${savecmd1}
    checkexec "unable to transfer result files from ${PWD}/ to ${resultdir}/" "succesfuly transferred result files from ${PWD}/ to ${resultdir}/"
  else
    ls ${dnchain}/${nfrad}* > /dev/null
    if [ ${?} == 0 ] ; then
      savecmd2="$savecmd ${dnchain}/${nfrad}* $resultdir/"
      echo "# $savecmd2"
      $savecmd2
      checkexec "unable to save result files from $dnchain to $resultdir/" "succesfuly transferred result files from $dnchain to $resultdir/"
    else
      echo "ERROR: unable to find the result files"
      exit 1
    fi
  fi
  if [[ "$worklocal" == "yes" ]] ; then
    # remove local copies of input/output files
    rm -f ./${nfrad}*
  fi
  
  echo ""
  echo "# # # # #"
  echo " # # # #"
done
