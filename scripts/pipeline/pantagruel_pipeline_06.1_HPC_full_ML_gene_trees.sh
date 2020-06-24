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
export defncpus='auto-configure'
export defmem='auto-configure'
export defwth='auto-configure'
export defhpctype='PBS'
export defchunksize=1000
export withpython='true'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

#############################################################
## 06.1 Full ML gene trees on HPC
#############################################################

mkdir -p ${ptglogs}/raxml/gene_trees/ ${mlgenetrees}/

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
sortminsizes=($(for famset in `ls ${genetrees}/*cdsfams_* | grep -v 'aln_list'` ; do echo $famset | sed -e 's/.*cdsfams_minsize\([0-9]\+.*\)/\1/' ; done | sort -n))
allcdsfam2phylo=($(for famminsize in ${sortminsizes[@]} ; do ls ${genetrees}/*${famminsize} ; done))
  
allmems=(4 8 32 64)
allwalltimes=(12 24 48 72)
allncpus=(4 4 8 16)

m=${mem}
w=${wth}
c=${ncpus}

for i in ${!allcdsfam2phylo[@]} ; do
  # iterate over lists of family classed by sizes and find matching setting
  cdsfam2phylo=${allcdsfam2phylo[$i]}
  [ ${i} -gt 3 ] && i=3 # max for settings
  # if their value was left to their default 'auto-configure', these variables are defined based on the size of the gene family
  [ "${m}" == 'auto-configure' ] && mem=${allmems[$i]}
  [ "${w}" == 'auto-configure' ] && wth=${allwalltimes[$i]}
  [ "${c}" == 'auto-configure' ] && ncpus=${allncpus[$i]}
  echo "cdsfam2phylo=${cdsfam2phylo} ; mem_per_core=${mem}gb ; walltime=${wth}:00:00 ; ncpus=${ncpus}"
  tasklist=${cdsfam2phylo}_aln_list
  if [[ "${resumetask}" == "true" && -s ${tasklist} ]] ; then
    echo "will re-use previously established list of gene families for which to compute trees:"
    ls -l ${tasklist}
  else
    rm -f ${tasklist}
    for fam in `cut -f1 $cdsfam2phylo` ; do ls ${cdsalifastacodedir}/${fam}.codes.aln >> $tasklist ; done
  fi
  if [ "$(wc -l $cdsfam2phylo | cut -f1 -d' ')" -lt "$(wc -l $tasklist | cut -f1 -d' ')" ] ; then 
    >&2 echo "ERROR $(dateprompt): missing gene family alignments; please fix the list '$tasklist' or the content of folder '$alifastacodedir/' before continuing."
    exit 1
  fi

  rm -f ${tasklist}*
  for fam in $(cut -f1 ${cdsfam2phylo}) ; do
    aln=${cdsalifastacodedir}/${fam}.codes.aln
    if [[ "${resumetask}" == "true" ]] ; then
      mainres=${mlgenetrees}/${mainresulttag}/RAxML_${mainresulttag}.${fam}.codes
      skiptag=${mlgenetrees}/bulk/${fam}.codes.smallreducedali
      if [[ ! -s ${mainres} && ! -s ${skiptag} ]] ; then
        ls ${aln} >> ${tasklist}_resume
      fi
    fi
    ls ${aln} >> ${tasklist}
  done
  if [[ "${resumetask}" == "true" ]] ; then
    if [[ -s ${tasklist}_resume ]] ; then
      echo "Resume task 6: $(wc -l ${tasklist}_resume | cut -d' ' -f1) ML trees left to infer"
    else
      echo "Resume task 6: all ML tree built; skip ML tree building"
    fi
    tasklist=${tasklist}_resume
  fi

  qsubvars="tasklist=${tasklist},outputdir=${mlgenetrees},reducedaln=true,nbthreads=${ncpus}"
  [ ! -z "${raxmlbin}" ] && qsubvars="${qsubvars},raxmlbin=${raxmlbin}"
  [ ! -z "${fwdenv}" ] && qsubvars="${qsubvars},${fwdenv}"
  [ ! -z "${modulefile}" ] && qsubvars="${qsubvars},modulefile=${modulefile}"
  
  [ -s ${tasklist} ] && Njob=`wc -l ${tasklist} | cut -f1 -d' '` || Njob=0
  [ ! -z "${topindex}" ] && [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
  
  # accomodate with possible upper limit on number of tasks in an array job
  jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
  for jobrange in ${jobranges[@]} ; do
    echo "jobrange=$jobrange"
    case "$hpctype" in
      'PBS') 
        qsub -J $jobrange -l walltime=${wth}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) \
	    -o ${ptglogs}/raxml/gene_trees -j oe -v "$qsubvars" ${ptgscripts}/raxml_array_PBS.qsub
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
	    nflog="${ptglogs}/raxml/gene_trees/raxml_gene_trees_$(basename $cdsfam2phylo).%J.%I.o"
	    bsub -J "raxml_gene_trees_$(basename $cdsfam2phylo)[$jobrange]" -q ${bqueue} -R "select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]" \
	    -n${ncpus} -M${memmb} -o ${nflog} -e ${nflog} -env "$qsubvars" ${ptgscripts}/raxml_array_LSF.bsub
	    ;;
	  *)
	    echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	    exit 1;;
    esac
  done
done

