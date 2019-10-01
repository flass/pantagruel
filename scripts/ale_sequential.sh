#!/bin/bash

tasklist=${1}
resultdir=${2}
spetree=${3}
nrecs=${4}
alealgo=${5}
maxpcmem=${6}

echo "ale_sequential.sh call was: ${@}"
echo "with variables set as: tasklist=${tasklist}, resultdir=${resultdir}, spetree=${spetree}, nrecs=${nrecs}, alealgo=${alealgo}, maxpcmem=${maxpcmem}"

usage() {
	echo "Usage: ale_sequential.sh tasklist resultdir spetree [nrecs; default=1000] [alealgo; default='ALEml_undated'] [max %mem; default=90]"
}

nbtaxafromtree () {
python << EOF
with open('${1}', 'r') as fchain:
  chainone = fchain.readline()
  print 'ntaxa:', chainone.count('(') + 2

EOF
}

cd ${ptgtmp}

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

# alealgo
echo "alealgo:"
if [ -z "$alealgo" ] ; then
  echo -n "Default: "
  alealgo='ALEml'
fi
echo "will use $alealgo algorithm for reconciliation estimation"
# relburninfrac
echo "relburninfrac:"
if [ -z "$relburninfrac" ] ; then
  echo -n "Default: "
  relburninfrac=0.25
fi
echo "will discard $relburninfrac fraction of the tree chain as burn-in"

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
if [ -z "$maxpcmem" ] ; then
  maxpcmem=90
fi

# worklocal
# worklocal='yes' indicates that there will be a specific machine where the 'compute work'
# is to be done 'locally' (e.g. if using a HPC system, a worker node) 
# that is different from the machine where the master/submission script is executed (e.g. login node),
# implying file traficking between them at the begin and end of the job.
echo "# worklocal:"
if [ -z "$worklocal" ] ; then
  echo "(Use default)"
  worklocal="no"
else
  if [[ "$worklocal"=="n" || "$worklocal"=="false" ]] ; then
    worklocal="no"
  elif [[ "$worklocal"=="y" || "$worklocal"=="true" ]] ; then
    worklocal="yes"
  fi
fi
echo "will work (read/write) locally: ${worklocal}"
echo ""
echo "# # # #"
echo ""



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

  # resume from run with already estimated parameters,
  # to perform further reconciliation sampling
  # (also allows to sample from defined parameter set)
  if [ "${#DTLrates[@]}" -eq 3 ] ; then
    echo -e "will perform analysis with set DTL rate parameters:\n${DTLrates[@]}"
  elif [ ! -z "$resumealefromtag" ] ; then
    estparam=($(ls ${resultdir}/${bnchain}.ale.*ml_rec${resumealefromtag}))
    if [ ! -z "${estparam}" ] ; then
    DTLrates=($(grep -A 1 "rate of" ${estparam} | grep 'ML' | awk '{print $2,$3,$4}'))
    if [ "${#DTLrates[@]}" -eq 3 ] ; then
      echo -e "will resume analysis from previously estimated DTL rate parameters:\n${DTLrates[@]} \nas found in file:\n'$estparam'"
      prevcomputetime=$(cat ${resultdir}/${nfrad}.ale.computetime${resumealefromtag} | cut -f3)
      if [ ! -z $prevcomputetime ] ; then
      echo -e "will add previous computation time spent estimating parameters found in file:\n'${dnchain}/${nfrad}.ale.computetime${resumealefromtag}'\nto new record:\n'./${nfrad}.ale.computetime'"
      fi
    fi
    fi
  fi
  echo ""

  if [[ -e ${nfchain}.ale ]] ; then
    if [[ "$worklocal"=="yes" ]] ; then
     # copy input files locally
     rsync -az ${nfchain}.ale ./
    fi
    echo "use pre-existing ALE index file:"
    ls ${nfchain}.ale
  elif [[ -e ${resultdir}/${bnchain}.ale ]] ; then
    if [[ "$worklocal"=="yes" ]] ; then
     # copy input files locally
     rsync -az ${resultdir}/${bnchain}.ale ./
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
  date

  # start timing in seconds
  SECONDS=0
  aleexe="${alebin}${alealgo}"
  # run ALE reconciliation 
  if [ "$alealgo" == 'ALEml' ] ; then
    alecmd="${aleexe} ${stree} ${chain}.ale ${nrecs} _"
    if [ "${#DTLrates[@]}" -eq 3 ] ; then alecmd="${alecmd} ${DTLrates[@]}" ; fi
  elif [ "$alealgo" == 'ALEml_undated' ] ; then
    alecmd="${aleexe} ${stree} ${chain}.ale sample=${nrecs} separators=_"
    if [ "${#DTLrates[@]}" -eq 3 ] ; then alecmd="${alecmd} delta=${DTLrates[0]} tau=${DTLrates[0]} lambda=${DTLrates[0]}" ; fi
  else
    echo "ALE algorithm $alealgo not supported in this script, sorry; exit now"
    exit 2
  fi
  echo "# ${alecmd}"
  # run it in bg with a daemon checking
  ${alecmd} &
  alepid=$!
  while [ ! -z $(ps -q ${alepid} -o comm=) ] ; do
    ## check memory use is not going off the charts
    ALEMEM=$(pmap ${alepid} | tail -n1 | awk '{print $NF}')
    echo "$nfrad\t$alealgo\t$ALEMEM\tkB" > ${nfrad}.ale.memusage
    ALETIME=$SECONDS
    echo -e "$nfrad\t$alealgo\t$ALETIME\ts" > ${nfrad}.ale.computetime
    pcmem=$(ps -o pid,%mem | grep ${alepid} | awk '{print $NF}')
    if [ ${pcmem%.*} -ge ${maxpcmem} ] ; then
      # stop immediately
      echo "Memory use is beyond ${maxpcmem}% the server capacity; stop the job now"
      kill -9 ${alepid}
    fi
    $savecmd ./${nfrad}.ale.* $resultdir/
    sleep 60s
  done
  
  echo ""
  echo "# # # #"

  ALETIME=$SECONDS
  if [ ! -z "$prevcomputetime" ] ; then ALETIME=$(( $ALETIME + $prevcomputetime )) ; fi
  echo -e "$nfrad\t$alealgo\t$ALETIME" > $nfrad.ale.computetime
  echo "reconciliation estimation took" $(date -u -d @${ALETIME} +"%Hh%Mm%Ss") "total time"
  if [ ! -z "$prevcomputetime" ] ; then echo "(including ${prevcomputetime} in previous run)" ; fi

  echo "# ls ./${nfrad}*"
  ls ./${nfrad}*

  # save files
  ls ./${nfrad}*.ale.* > /dev/null
  if [ $? == 0 ] ; then
    savecmd1="$savecmd ./${nfrad}*.ale* $resultdir/"
    echo "# $savecmd1"
    $savecmd1
    checkexec "unable to transfer result files from $PWD/ to $resultdir/" "succesfuly transferred result files from $PWD/ to $resultdir/"
  else
    ls ${dnchain}/${nfrad}*.ale.* > /dev/null
    if [ $? == 0 ] ; then
      savecmd2="$savecmd ${dnchain}/${nfrad}*.ale.* $resultdir/"
      echo "# $savecmd2"
      $savecmd2
      checkexec "unable to save result files from $dnchain to $resultdir/" "succesfuly transferred result files from $dnchain to $resultdir/"
    else
      echo "ERROR: unable to find the result files"
      exit 1
    fi
  fi
  if [[ "$worklocal"=="yes" ]] ; then
    # remove local copies of input/output files
    rm -f ./${nfrad}*
  fi
  
  echo ""
  echo "# # # # #"
  echo " # # # #"
done
