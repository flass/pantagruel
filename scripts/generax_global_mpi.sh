#!/bin/bash
### verify and parse key variable definition
# generaxfamfi
if [ -z "${generaxfamfi}" ] ; then
  generaxfamfi="${1}"
  if [ -z "${generaxfamfi}" ] ; then
    echo "missing mandatory argument: generax families file path (to pass either through env var \$generaxfamfi or as positional argument \$1); exit now"
	exit 1
  fi
fi
if [ -z "${outrecdir}" ] ; then
  outrecdir="${2}"
  if [ -z "${outrecdir}" ] ; then
    echo "missing mandatory argument: generax output folder path (to pass either through env var \$outrecdir or as positional argument \$2); exit now"
	exit 1
  fi
fi
mkdir -p ${outrecdir}/
# spetree
echo "spetree:"
if [ -z "${spetree}" ] ; then
  spetree=${speciestree}_clade_defs.nwk
  echo "default value for species tree: ${spetree}"
fi
if [ ! -s "${spetree}" ] ; then
  echo "Error: species tree file '$spetree' does not exist/is empty; exit now"
  exit 2
fi
if [ "${GeneRaxalgo}" == 'reconciliation-samples' ] ; then
  generaxopt="${generaxopt} --reconciliation-samples ${recsamplesize}"
else
  generaxopt="${generaxopt} --reconcile"
fi
# generaxcommonopt
if [ -z "${generaxcommonopt}" ] ; then
  generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR" # env var is a pipeline default
fi
# ncpus
if [ -z "${ncpus}" ] ; then
  ncpus=1
fi

if [ ${ncpus} -gt 1 ] ; then
  grxcmd="mpirun -np ${ncpus} generax"
else
  grxcmd="generax"
fi

${grxcmd} ${generaxcommonopt} -s ${spetree} -f ${generaxfamfi} -p ${outrecdir} ${generaxopt}
