#!/bin/bash

### verify and parse key variable definition
# generaxfamfi
if [ -z ${generaxfamfi} ] ; then
  generaxfamfi="${1}"
  if [ -z ${generaxfamfi} ] ; then
    echo "missing mandatory argument: generax family file path (to pass either through env var \$generaxfamfi or as positional argument \$1); exit now"
	exit 1
  fi
fi
# outrecdir
echo "outrecdir:"
if [ -z "${outrecdir}" ] ; then
  echo "ERROR: need to define variable outrecdir ; exit now"
  exit 2
else
  ls ${outrecdir} -d
  if [ ${?} != 0 ] ; then
    echo "directory '${outrecdir}' is missing ; create it now"
    mkdir -p ${outrecdir}/
    if [ ${?} != 0 ] ; then
      echo "could not create directory '${outrecdir}' ; exit now"
      exit 2
    fi
  fi
fi
# spetree
echo "spetree:"
if [ -z "${spetree}" ] ; then
  echo "ERROR: need to define variable spetree ; exit now"
  exit 2
fi
# generaxcommonopt
if [ -z "${generaxcommonopt}" ] ; then
  generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR" # env var is a pipeline default
fi
# GeneRaxalgo
if [ ${GeneRaxalgo} == 'reconciliation-samples' ] ; then
  generaxopt="${generaxopt} --reconciliation-samples ${recsamplesize}"
fi
# ncpus
if [ ! -z "${ncpus}" ] ; then
  ncpus=1
fi

######

nfrad1=$(basename ${generaxfamfi})
dngrff=$(dirname ${generaxfamfi})
nfext=${nfrad1##*.}
nfrad2=${nfrad1%.*}
echo ${nfrad2}

ls ${spetree}
if [ ${?} != 0 ] ; then
  echo "look for ${spetree} species tree file in ${dnchain}/ folder"
  ls ${dngrff}/${nfrad}*${spetree}*
  if [ ${?} != 0 ] ; then
      echo "ERROR: file '${spetree}' is missing/empty ; exit now"
      exit 2
  else
    echo "found it!" 
    lnfstree=(`ls ${dngrff}/${nfrad}*${spetree}* 2> /dev/null/`)
    nfstree=${lnfstree[0]}
    echo "will use nfstree=${nfstree}"
  fi
else
  nfstree=${spetree}
fi

mkdir -p ${outrecdir}

if [ ${ncpus} -gt 1 ] ; then
  grxcmd="mpirun -np ${ncpus} generax"
else
  grxcmd="generax"
fi
mpirun -np ${ncpus} generax ${generaxcommonopt} -s ${speciestree}_clade_defs.nwk -f ${generaxfamfi} -p ${outrecdir}/${nfrad2} ${generaxopt}