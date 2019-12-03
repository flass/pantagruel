#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 6 Nov 2019

# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=8
export defmem=32
export defwth=24
export defhpctype='PBS'
export defchunksize=3000
export withpython='true'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

#############################################################
## 06.2 Gene tree collapsing on HPC
#############################################################

# first collate the list of families for which no gene tree was computed 
# due to containing too few non-identical sequences 
# and will be discarded for the remamining analyses
smallfams=${genetrees}/reduced_alignment_is_too_small_fams
touch ${smallfams}
cat ${genetrees}/bulk/*.smallreducedali > ${smallfams}

echo "$(wc -l ${smallfams}) gene families were found to contain to few non-identical CDSs for a gene tree to be computed"

#### OPTION: edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # no collapsing, just convert the alignments from fasta to nexus
  for aln in `ls ${alifastacodedir}` ; do
    convalign -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  
else
 if [[ "${chaintype}" == 'collapsed' ]] ; then
  if [ -z ${collapsecolid} ] ; then
    collapsecolid=1
  fi
  if [[ "$collapseCladeParams" != 'default' ]] ; then
    eval "$collapseCladeParams"
    # e.g.:  eval 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'
  fi
  
  ## detect clades to be collapsed in gene trees
  collapsecriteriondef="--clade_stem_conds=[('$criterion','>=',$cladesupp)] --within_clade_conds=[('$withinfun','$criterion','<=',$subcladesupp,-1),('max','$criterion','<',$cladesupp,-1)]"
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  echo "${collapsecriteriondef}" > ${colalinexuscodedir}/${collapsecond}.collapse_criterion_def
  mlgenetreelist=${mlgenetrees%*s}_list
  ${ptgscripts}/lsfullpath.py "${mlgenetrees}/${mainresulttag}/*" > ${mlgenetreelist}
  
  # accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
  Njob=`wc -l ${mlgenetreelist} | cut -f1 -d' '`
  [ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
  
  jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))

  # In the job submission commands below, some lines are specific to the HPC system 
  # and environmnt on which the script was developped:
  #   module load anaconda3/personal
  #   source activate env_python2
  # In other environments, other methods may be used to access the required Python packages.
  # To emulate this on other systems, it is best to use anaconda to build your own environment
  # where you make sure you are using python 2.7 and have all required packages installed with it.
  # You can do so using the following command:
  # conda create -n env_python2 python=2.7 python-igraph biopython bcbio-gff scipy

  for jobrange in ${jobranges[@]} ; do
    beg=`echo $jobrange | cut -d'-' -f1`
    tail -n +${beg} ${mlgenetreelist} | head -n ${chunksize} > ${mlgenetreelist}_${jobrange}
	
    qsubvars="mlgenetreelist=${mlgenetreelist}_${jobrange},ptgscripts,cdsalifastacodedir,ncpus,colalinexuscodedir,collapsecond,mlgenetrees"
    if [ ! -z "${fwdenv}" ] ; then
	  qsubvars="${qsubvars},${fwdenv}"
	fi
	
	case "${hpctype}" in
      'PBS')
         subcmd="qsub -N 'mark_unresolved_clades' -l select=1:ncpus=${ncpus}:mem=${mem}gb,walltime=${wth}:00:00 \
          -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe \
		  -v \"${qsubvars}\" ${ptgscripts}/mark_unresolved_clades.qsub"
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
	    nflog="${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.o"
	    subcmd="bsub -J 'mark_unresolved_clades' -q ${bqueue} \
	     -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" \
         -n${ncpus} -M${memmb} -env \"$qsubvars\" \
         -o ${nflog} -e ${nflog} ${ptgscripts}/mark_unresolved_clades.bsub"
	    ;;
	  *)
	    echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	    exit 1;;
    esac
	echo "${subcmd}"
	eval "${subcmd}"
  done
  export collapsecoldate=$(date +%Y-%m-%d)
  echo -e "${collapsecolid}\t${collapsecoldate}\t${collapsecriteriondef}" > ${genetrees}/collapsecol
  
 else
  echo "Error: incorrect value for variable chaintype: '${chaintype}'"
  exit 1
 fi
fi
#### end OPTION
