#!/bin/bash

tasklist="${1}"
export outputdir="${2}"
export mbmcmcopt="${3}"

if [ -z $q ] ; then
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

### e.g. tasklist=/path/list  ; outputdir=/path/dir ; mbmcmcopt="params string"

# test key variables have been provided
if [ -z "${tasklist}" ] ; then
  echo "!!! ERROR : mandatory variable \${tasklist} not declared ; exit now"
  exit 1
fi
if [ -z "${outputdir}" ] ; then
  echo "!!! ERROR : mandatory variable \${outputdir} not declared ; exit now"
  exit 1
fi
# $mbmcmcopt and $mbmcmcpopt can be empty, which will make MrBayes run defaut parameters



mkdir -p ${outputdir}/
cd ${outputdir}/
if [ ${?} != 0 ] ; then 
  echo "!!! ERROR : unable to access output directory '${outputdir}/' ; exit now"
  exit 1
fi
if [ ! -e "${tasklist}" ] ; then 
  echo "!!! ERROR : unable to access task list file '${tasklist}' ; exit now"
  exit 1
fi



runmb(){

nfaln=$1
nfrad1=$(basename ${nfaln})
nfrad2=${nfrad1%.*}
echo ${nfrad2}
echo ""

pref=$(echo ${nfrad2} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}")
mkdir -p ${outputdir}/${pref}/
cd ${outputdir}/${pref}/


echo "current directory (output directory) is $HOSTNAME:$PWD"

echo "mbmcmcopt=$mbmcmcopt"
#
#mbresume=$(python2.7 << EOF
#mbmcmcopts="$mbmcmcopt".lower().split()
#mbresume = 'no'
#for opteqval in mbmcmcopts:
#  opt, val = opteqval.split('=')
#  if opt.strip(' "')=='append': mbresume = val.strip(' "')
#print mbresume
#EOF
#)
if [ "${mbresume}" == 'yes' ] ; then
  if [ "$PWD" != "$(realpath ${outputdir}/${pref})" ] ; then
    # import files from previous interupted analysis
    echo "rsync -avz ${outputdir}/${pref}/*${nfrad2}* ./"
    rsync -avz ${outputdir}/${pref}/*${nfrad2}* ./
	echo "recovered intermediary files from previous run with exit status $?"
	[ -s ./*${nfrad2}*ckp ] && mbmcmcopt="$mbmcmcopt append=yes"
  fi
else
  if [ ! -z "$(ls ./${nfrad2}.mb* 2> /dev/null)" ] ; then
    echo "there are already files matching pattern './${nfrad2}.mb*' ; I do not risk overwritting them: exit now"
    exit 1
  fi
fi
echo "rsync -avz ${nfaln} ./"
rsync -avz ${nfaln} ./
echo "copied input files with exit status $?"
echo "ls ./*${nfrad2}*"
ls ./*${nfrad2}*
echo ""

# give a glimpse of data complexity
echo "data matrix:"
\grep 'dimensions' ${nfrad1}
echo ""

# set MrBayes parameters
echo "set autoclose=yes nowarn=yes" > ${nfrad2}.mbparam.txt
echo "execute ${nfrad1}" >> ${nfrad2}.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
  echo "lset nst=6 rates=invgamma ploidy=haploid" >> ${nfrad2}.mbparam.txt
else
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 due to a bug not recognizing ploidy option
  echo "lset nst=6 rates=invgamma" >> ${nfrad2}.mbparam.txt
fi
echo "showmodel" >> ${nfrad2}.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 as the checkpoint option is new to v3.2
  echo "mcmcp $mbmcmcopt file=${nfrad2}.mb checkpoint=yes checkfreq=100000" >> ${nfrad2}.mbparam.txt
else
  echo "mcmcp $mbmcmcopt file=${nfrad2}.mb" >> ${nfrad2}.mbparam.txt
fi
echo "mcmc" >> ${nfrad2}.mbparam.txt
if [ "${mbversion%.*}" == '3.2' ] ; then
# following line WORKED with MrBayes v3.2.2, but DID NOT WORK with v3.1.2 as the conformat option is new to v3.2
  echo "sumt conformat=simple contype=allcompat minpartfreq=0.001" >> ${nfrad2}.mbparam.txt
else
  echo "sumt minpartfreq=0.001 contype=allcompat" >> ${nfrad2}.mbparam.txt
fi
echo "sump" >> ${nfrad2}.mbparam.txt
echo "quit" >> ${nfrad2}.mbparam.txt

cat ${nfrad2}.mbparam.txt
echo ""


# run the progrram
echo "running MrBayes:"
echo "mb < ${nfrad2}.mbparam.txt"
mb < ${nfrad2}.mbparam.txt > ${nfrad2}.mb.log


echo "output of MrBayes phylogenetic reconstruction is :"
echo "ls ./*${nfrad2}*"
ls ./*${nfrad2}*
echo ""


echo ""
echo "- - - - - - - - - - - - - - - - - - "
echo ""
}

export -f runmb 
# choose to run each tree on only one process, and to parallelize by the tasks (avoid relying on MPI interface)
parallel runmb :::: ${tasklist}
