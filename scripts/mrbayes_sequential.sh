#!/bin/bash

tasklist="${1}"
outputdir="${2}"
export mbmcmcopt="${3}"

if [ -z $mbversion ] ; then
# assumes use of most recent version
mbversion='3.2.6'
fi

## I/O file path (lists)
# tasklist
# outputdir
## the MCMCMC parameters, including number of parallel runs and number of mixed chains; string to be passed to mcmcp command, e.g. : "Nruns=2 Nchains=6 Ngen=2000000"
# mbmcmpcopt
## the MCMCMC run options, notably the flag for resuming from a checkpoint
# mbmcmcopt

### e.g. qsub -v tasklist=/path/list outputdir=/path/dir mbmcmcopt="params string"

# test key variables have been provided
if [ -z "$tasklist" ] ; then
  echo "!!! ERROR : mandatory variable \$tasklist not declared ; exit now"
  exit 1
fi
if [ -z "$outputdir" ] ; then
  echo "!!! ERROR : mandatory variable \$outputdir not declared ; exit now"
  exit 1
fi
# $mbmcmcopt and $mbmcmcpopt can be empty, which will make MrBayes run defaut parameters


mkdir -p $outputdir/
cd $outputdir/
if [ $? != 0 ] ; then 
  echo "!!! ERROR : unable to access output directory '$outputdir/' ; exit now"
  exit 1
fi
if [ ! -e "$tasklist" ] ; then 
  echo "!!! ERROR : unable to access task list file '$tasklist' ; exit now"
  exit 1
fi



runmb(){

nfaln=$1
nfrad1=$(basename $nfaln)
nfrad2=${nfrad1%.*}
echo $nfrad2
echo ""

echo "current directory (output directory) is $HOSTNAME:$PWD"

#~ # remove traces from potential previous chains 
#~ rm ./*$nfrad2*
#~ echo "removed pre-existing files with name containing '$nfrad2' with exit status $?"
#~ echo ""
# directory is supposed to be clean anyway

echo 'mbmcmcopt="$mbmcmcopt"'

mbresume=$(python << EOF
mbmcmcopts="$mbmcmcopt".lower().split()
mbresume = 'no'
for opteqval in mbmcmcopts:
  opt, val = opteqval.split('=')
  if opt.strip(' "')=='append': mbresume = val.strip(' "')
print mbresume
EOF
)
if [ "$mbresume" == 'yes' ] ; then
  # import files from previous interupted analysis
  echo "rsync -avz $outputdir/*${nfrad2}* ./"
  rsync -avz $outputdir/*${nfrad2}* ./
else
  echo "rsync -avz $nfaln ./"
  rsync -avz $nfaln ./
fi
echo "copied input files with exit status $?"
echo "ls ./"
ls ./
echo ""

# give a glimpse of data complexity
echo "data matrix:"
\grep 'dimensions' $nfrad1
echo ""

# set MrBayes parameters
echo "set autoclose=yes nowarn=yes" > $nfrad2.mbparam.txt
echo "execute $nfrad1" >> $nfrad2.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
  echo "lset nst=6 rates=invgamma ploidy=haploid" >> $nfrad2.mbparam.txt
else
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 due to a bug not recognizing ploidy option
  echo "lset nst=6 rates=invgamma" >> $nfrad2.mbparam.txt
fi
echo "showmodel" >> $nfrad2.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 as the checkpoint option is new to v3.2
  echo "mcmcp $mbmcmcopt file=$nfrad2.mb.nex checkpoint=yes checkfreq=100000" >> $nfrad2.mbparam.txt
else
  echo "mcmcp $mbmcmcopt file=$nfrad2.mb.nex" >> $nfrad2.mbparam.txt
fi
echo "mcmc" >> $nfrad2.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 as the conformat option is new to v3.2
  echo "sumt conformat=simple contype=allcompat minpartfreq=0.001" >> $nfrad2.mbparam.txt
else
  echo "sumt minpartfreq=0.001 contype=allcompat" >> $nfrad2.mbparam.txt
fi
echo "sump" >> $nfrad2.mbparam.txt
echo "quit" >> $nfrad2.mbparam.txt

cat $nfrad2.mbparam.txt
echo ""


# run the progrram in background
echo "running MrBayes:"
#~ echo "mpirun -np $nbcores mb < $nfrad2.mbparam.txt &"
#~ mpirun -np $nbcores mb < $nfrad2.mbparam.txt &
#~ echo "mpiexec mb < $nfrad2.mbparam.txt"
#~ mpiexec mb < $nfrad2.mbparam.txt
echo "mb < $nfrad2.mbparam.txt"
mb < $nfrad2.mbparam.txt


echo "output of MrBayes phylogenetic reconstruction is :"
echo "ls ./*$nfrad2*"
ls ./*$nfrad2*
echo ""


echo ""
echo "- - - - - - - - - - - - - - - - - - "
echo ""
}

export -f runmb 
# choose to run each tree on only one process, and to parallelize by the tasks (avoid relying on MPI interface)
parallel runmb ::: `cat $tasklist`
