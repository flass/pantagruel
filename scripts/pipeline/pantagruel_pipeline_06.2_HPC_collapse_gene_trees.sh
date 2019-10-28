#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [ncpus=8] [mem=32gb] [walltimehours=4] [hpc_type:{PBS(default)|LSF}]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  ncpus="$2"
else
  ncpus=8
fi
if [ ! -z "$3" ] ; then
  mem="$3"
else
  mem=32
fi
if [ ! -z "$4" ] ; then
  wth="$4"
else
  wth=4
fi
if [ -z "$5" ] ; then 
  hpctype='PBS'
else
  hpctype="$5"
fi

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
  #~ export nexusaln4chains=${colalinexuscodedir}
  #~ export mboutputdir=${bayesgenetrees}
  
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
  #~ export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
  #~ export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
  collapsecriteriondef="--clade_stem_conds=[('$criterion','>=',$cladesupp)] --within_clade_conds=[('$withinfun','$criterion','<=',$subcladesupp,-1),('max','$criterion','<',$cladesupp,-1)]"
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  echo "${collapsecriteriondef}" > ${colalinexuscodedir}/${collapsecond}.collapse_criterion_def
  mlgenetreelist=${mlgenetrees%*s}_list
  ${ptgscripts}/lsfullpath.py "${mlgenetrees}/${mainresulttag}/*" > ${mlgenetreelist}
  
  # accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
  Njob=`wc -l ${mlgenetreelist} | cut -f1 -d' '`
  chunksize=3000
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
    qsubvars="ptgscripts=${ptgscripts},mlgenetreelist=${mlgenetreelist}_${jobrange},cdsalifastacodedir=${cdsalifastacodedir},ncpus=${ncpus},colalinexuscodedir=${colalinexuscodedir},collapsecond=${collapsecond},didseq=${mlgenetrees}/identical_sequences"
    case "$hpctype" in
      'PBS')
         qsub -N 'mark_unresolved_clades' -l select=1:ncpus=${ncpus}:mem=${mem}gb,walltime=${wth}:00:00 \
          -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe \
		  -v "${qsubvars}" ${ptgscripts}/mark_unresolved_clades.qsub
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
	    bsub -J 'mark_unresolved_clades' -q ${bqueue} \
	     -R "select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]" \
         -n${ncpus} -M${memmb} -env "$qsubvars" \
         -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.o \
         -e ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.e \
         ${ptgscripts}/mark_unresolved_clades.bsub
	    subcmd="bsub -J 'mark_unresolved_clades' -q ${bqueue} \
	     -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" \
         -n${ncpus} -M${memmb} -env \"$qsubvars\" \
         -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.o \
         -e ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.e \
         ${ptgscripts}/mark_unresolved_clades.bsub"
		 echo ${subcmd}
	    ;;
	  *)
	    echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	    exit 1;;
    esac
  done
  export collapsecoldate=$(date +%Y-%m-%d)
  echo -e "${collapsecolid}\t${collapsecoldate}" > ${genetrees}/collapsecol
  #~ export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
  #~ export mboutputdir=${bayesgenetrees}/${collapsecond}
  
 else
  echo "Error: incorrect value for variable chaintype: '${chaintype}'"
  exit 1
 fi
fi
#### end OPTION
