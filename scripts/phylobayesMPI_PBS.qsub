#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=2:ncpus=32:mem=64gb
#PBS -N phylobayesMPI
#PBS -o /home/flassall/logs/raxml/gene_trees
#PBS -j oe

# start timing
SECONDS=0

## load modules
module load intel-suite
module load gsl

echo "passed variables:"
echo "  nfaln=${nfaln}"
echo "  outputdir=${outputdir}"
echo "  chainname=${chainname}"
echo "  pbbin=${pbbin}"
echo "  coreperchain=${coreperchain}"
echo "  nbchain=${nbchain}"
echo "  pboptions=${pboptions}"
echo "  pbresume=${pbresume}"
echo "  timelimit=${timelimit}"
echo "  syncfreq=${syncfreq}"

#~ ## user-installed openmpi
#~ module load gcc/4.9.1
#~ # code to install the software = just indications to reproduce -- VERY LONG!
#~ # wget --no-check-certificate https://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.8.tar.gz
#~ # tar -xvzf openmpi-1.8.8.tar.gz
#~ # mkdir -p ${HOME}/software/openmpi-1.8.8-build ${HOME}/libraries/openmpi-1.8.8
#~ # cd ${HOME}/software/openmpi-1.8.8-build
#~ # ../openmpi-1.8.8/configure --prefix=${HOME}/libraries/openmpi-1.8.8/
#~ # make all install
#~ # load libs
#~ export LD_LIBRARY_PATH=${HOME}/libraries/openmpi-1.8.8/lib:${LD_LIBRARY_PATH}
#~ export PATH=${HOME}/libraries/openmpi-1.8.8/bin:${PATH}
#~ export CPATH=${HOME}/libraries/openmpi-1.8.8/include:${CPATH}
#~ export CPPPATH=${HOME}/libraries/openmpi-1.8.8/include:${CPPPATH}
#~ export C_INCLUDE_PATH=${HOME}/libraries/openmpi-1.8.8/include:${C_INCLUDE_PATH}
#~ export INCLUDE_PATH=${HOME}/libraries/openmpi-1.8.8/include:${C_INCLUDE_PATH}
#~ export MANPATH=${HOME}/libraries/openmpi-1.8.8/share/man:${MANPATH}

if [ -z $mpibin ] ; then
  mpibin='mpiexec'
  #~ mpibin='mpirun'
fi
if [ "$mpibin" == 'none' ] ; then
  echo "WILL NOT USE MPI"
else
  module load mpi
  module list
  echo "using $mpibin:"
  which $mpibin
  if [ $? != 0 ] ; then
    echo "!!! ERROR : unable to find $mpibin ; exit now"
    exit 1
  fi
fi

# test user-provided variables are correctly set
## mandatory variables
if [ -z $nfaln ] ; then
  echo "!!! ERROR : mandatory variable \$nfaln not declared ; exit now"
  exit 1
fi
if [ ! -e $nfaln ] ; then
  echo "!!! ERROR : cannot access mandatory alignment file '$nfaln' ; exit now"
  exit 1
fi
nfrad1=$(basename $nfaln)
nfrad2=${nfrad1%.*}
echo $nfrad2
if [ -z $outputdir ] ; then
  echo "!!! ERROR : mandatory variable \$outputdir not declared ; exit now"
  exit 1
fi
mkdir -p $outputdir/
  if [ ! -d $outputdir/ ] ; then 
    echo "!!! ERROR : unable to access output directory '$outputdir/' ; exit now"
    exit 1
fi

# variables with default values
# define software executables
# pb
if [ "$mpibin" != 'none' ] ; then
  if [ -z "$pbbin" ] ; then
    pbbin='/home/flassall/software/pbmpi/data/pb_mpi'
    echo "set default pbbin='$pbbin'"
  fi
  echo "will use Phylobayes MPI binary:"
  ls $pbbin
else
  if [ -z "$pbbin" ] ; then
    pbbin='/home/flassall/software/phylobayes4.1c/data/pb'
    echo "set default pbbin='$pbbin'"
  fi
  echo "will use Phylobayes (sequential version) binary:"
  ls $pbbin
fi
if [ $? != 0 ] ; then
  echo "!!! ERROR : unable to access Phylobayes MPI binary ; exit now"
  exit 1
fi
if [ ! -x $pbbin ] ; then
  echo "!!! ERROR : Phylobayes MPI binary '$pbbin' is not executable ; exit now"
  exit 1
fi

# bpcomp
if [ -z "$bpcompbin" ] ; then
  bpcompbin=$(dirname $pbbin)/bpcomp
fi
if [ -x "$bpcompbin"  ] ; then
  echo "will assess convergence with:"
  ls $bpcompbin
else
  echo "no bpcomp executable provided; will NOT assess convergence."
  bpcompbin=''
fi

# define parallelism parameters
if [ "$mpibin" != 'none' ] ; then
  if [ -z "$nbchain" ] ; then
    nbchain=2
  fi
  case ${nbchain} in
    0|*[!0-9]*) echo "!!! ERROR, number of chains '\$nbchain' must be non-null positive integer" ; exit 1;;
    *) echo "will run $nbchain separate chains" ;;
  esac
  if [ -z "$coreperchain" ] ; then
    coreperchain=8
  fi
  case ${coreperchain} in
    0|1|*[!0-9]*) echo "!!! ERROR, number of core chains '\$coreperchain' must be an integer >=2" ; exit 1;;
    *) echo "will run each chain on $coreperchain cores" ;;
  esac
fi

# define software options
if [ -z "$pboptions" ] ; then
  pboptions="-cat -gtr -dgam 4"
  echo "set default pboptions='$pboptions'"
fi
# define running time limit
if [ -z "$timelimit" ] ; then
  timelimit=$(( 72 * 60 * 60 ))
  echo "set default timelimit=$timelimit"
fi
# sync frequency
if [ -z "$syncfreq" ] ; then
  syncfreq=30m
  echo "set default syncfreq='$syncfreq'"
fi
syncf=`python2.7 -c "t = '$syncfreq' ; dmul = {'h':3600, 'm':60, 's':1} ; ts = int(t[:-1])*dmul[t[-1]] ; print ts"`
if [ $? != 0 ] ; then
  echo "!!! ERROR : wrong format for sync frequency: $syncfreq ; expect X{h,m,s} with X any integer ; exit now"
  exit 1
fi
# chain name
if [ -z "$chainname" ] ; then
  pchain=${nfrad2}.pb.chain
  echo "set default pchain='$pchain'"
else
  pchain=${chainname}
fi

# if a dry run
if  [ "$dryrun" == 'yes' ] ; then
  echo "# WOULD import files with this command:"
  echo "rsync -avz $nfaln $TMPDIR/"
  if [ "$pbresume" == 'yes' ] ; then
    echo "# resuming previously begun analysis; WOULD import files with this command:"
    echo "rsync -avz $outputdir/${pchain}* $TMPDIR/"
    echo "# check presence of previous run files"
    echo "ls $outputdir/${pchain}*"
    ls $outputdir/${pchain}*
  fi
  ### end of dry run
  echo "end of dry run ; exit now"
  exit 0
fi

## copy data locally
cd $TMPDIR/
echo "current directory is $PWD"
rsync -avz $nfaln $TMPDIR/
if [ $? != 0 ] ; then
  echo "!!! ERROR : unable to copy input file $nfaln into $TMPDIR/ ; exit now"
  exit 1
else
  echo "copied input files $nfaln succesfully"
fi
echo "ls ./"
ls ./

if [ "$pbresume" == 'yes' ] ; then
 echo "# resuming previously begun analysis; import files:"
  echo "rsync -avz $outputdir/${pchain}* $TMPDIR/"
  rsync -avz $outputdir/${pchain}* $TMPDIR/
  echo "copied previously computed chain files with exit status $?"
  echo "ls ./"
  ls ./
  echo ""
fi

# set PB command
if [ "$mpibin" == 'none' ] ; then
  pbcmd="${pbbin}" 
else
  pbcmd="$mpibin -n ${coreperchain} ${pbbin}"
fi
if [ "$pbresume" != 'yes' ] ; then
  pbcmd="${pbcmd} -d ${nfrad1} ${pboptions}"
  echo -e "##### Start PB chains\n"
else
  echo -e "##### Resume previous PB chains\n"  
fi

# start chain(s)
pbpid=()
if [ $nbchain -eq 1 ] ; then
  echo -e "${pbcmd} ${pchain} &> ${pchain}.log &\n"
  ${pbcmd} ${pchain} &> ${pchain}.log &
  pbpid+=($!)
else
  for i in $(seq 1 ${nbchain}) ; do
    echo -e "${pbcmd} ${pchain}.${i} &> ${pchain}.${i}.log &\n"
    ${pbcmd} ${pchain}.${i} &> ${pchain}.${i}.log &
    pbpid+=($!)
  done
fi
pbpidorstr=$(echo ${pbpid[@]} | sed -e 's/ /\\|/g')

# set and run daemon
rsynccmd="rsync -avz $TMPDIR/${pchain}* $outputdir/"
echo "will run '$rsynccmd' every ${syncfreq} (${syncf}s)"
lastsync=0
# run until the near end of time
while [[ ${SECONDS} -lt $(( ${timelimit} - (2 * ${timelimit} / 100) )) ]] ; do
  sleep 60s ; echo -n "."
  lastsync=$(( ${lastsync} + 60 ))
  # syncing and convergence check
  if [[ ${lastsync} -ge ${syncf} ]] ; then
    top -b -n1 -u $USER -b >> $TMPDIR/${pchain}.toplog
    echo ""
    $rsynccmd >> $TMPDIR/${pchain}.rsynclog
    lastsync=0
    # convergence diagnostic for possible early end of execution
    if [[ ! -z $bpcompbin && $nbchain -gt 1 ]] ; then
      nburnin=$(( $(wc -l ${pchain}.1.treelist | cut -d' ' -f1) / 4 ))
      $bpcompbin -x $nburnin $(ls ${pchain}.*.chain | sed -e 's/.chain$//g')
      cat bpcomp.bpdiff | xargs
      maxdiff=$(grep 'maxdiff' bpcomp.bpdiff | awk '{print $NF}')
      if [ -z $maxdiff ] ; then maxdiff=1 ; fi
      maxdiff10=$(python2.7 -c "md=float($maxdiff) ; print int(md*10)")
      if [ $maxdiff10 -lt 1 ] ; then
        # maxdiff < 0.1
        echo "convergence of topologies between chains detected, will stop the runs now"
        for fbpcomp in `ls bpcomp.*` ; do mv $fbpcomp ${pchain}.${fbpcomp#bpcomp.} ; done
        break
      fi
    fi
  fi
  # check the programs are still running...
  if [ -z "$(ps | grep $pbpidorstr)" ] ; then
    echo ""
    echo "!!! ERROR : unexpected termination of all pb chains ; perform last rsync and exit"
    $rsynccmd
    exit 1
  fi
done
echo ""
# final account of activity, broken-down by thread
top -H -b -n1 -u $USER -b >> $TMPDIR/${pchain}.toplog
# request chains to stop
for i in $(seq 1 ${nbchain}) ; do
  echo 0 > $TMPDIR/${pchain}.${i}.run
done
# wait for completion of chains' last cycles
while [ ! -z "$(ps | grep $pbpidorstr)" ] ; do
  sleep 60s
done
$rsynccmd

## end of script



