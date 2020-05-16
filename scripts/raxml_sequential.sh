#!/bin/bash

tasklist="${1}"
outputdir="${2}"
model="${3}"
keepresulttag="${4}"
bootstrapalgo="${5}"
nbthreads="${6}"
reducedaln="${7}"

# mandatory variables
if [ -z "$tasklist" ] ; then
  echo "!!! ERROR : must provide a task list through 'tasklist' variable ; exit now"
  exit 1
fi
if [ -z "${outputdir}" ] ; then
  echo "!!! ERROR : must provide a output directory path through 'outputdir' variable ; exit now"
  exit 1
fi
# options
if [ -z "${model}" ] ; then
  model='GTRCAT'
fi
if [ -z "${keepresulttag}" ] ; then
  resulttags=('bipartitions' 'rootedTree' 'identical_sequences')
else
  resulttags=("${keepresulttag}")
fi
if [ -z "${bootstrapalgo}" ] ; then
  bootstrapalgo='x'
fi
if [ -z "${nbthreads}" ] ; then
  nbthreads=$(nproc)
fi
if [ -z "${reducedaln}" ] ; then
  reducedaln=false
fi

for resulttag in bulk ${resulttags[@]} ; do
  mkdir -p ${outputdir}/$resulttag/
  if [ ! -d ${outputdir}/$resulttag ] ; then 
    echo "!!! ERROR : unable to access output directory 'outputdir/${resulttag}/' ; exit now"
    exit 1
  fi
done
if [ ! -e "${tasklist}" ] ; then 
  echo "!!! ERROR : unable to access task list file '${tasklist}' ; exit now"
  exit 1
fi

raxlogs=${ptglogs}/raxml/gene_trees
mkdir -p ${raxlogs}

### begin loop over alignment tasks
for nfaln in `cat ${tasklist}` ; do

raxlog=${raxlogs}/${nfrad2}.raxml.log
echo -e "# this is Pantagruel script raxml_sequential.sh\ntask: '${nfaln}'; current directory:${PWD} ; date: $(date "+%Y-%m-%d_%H-%M-%S")" &> ${raxlog}

diraln=$(dirname ${nfaln})
nfrad1=$(basename ${nfaln})
nfext=${nfrad1##*.}
nfrad2=${nfrad1%.*}
echo "# ouput will be named with file radical nfrad2='${nfrad2}'" &>> ${raxlog}


cd ${outputdir}/bulk/

# convert file format
if [ "${nfext}" == 'nex' ] ; then
  python2.7 -c "from Bio import AlignIO ; AlignIO.convert('${nfaln}', 'nexus', '${outputdir}/bulk/${nfrad2}.fasta', 'fasta')" &>> ${raxlog}
  if [ $? == 0 ] ; then
    echo "succesfully converted Nexus input file ${nfrad1} into FASTA format: ${nfrad2}.fasta" &>> ${raxlog}
    localn=${outputdir}/bulk/${nfrad2}.fasta
  else
    echo "failed conversion of input file ${nfrad1} into FASTA format; we'll see if it goes through..." &>> ${raxlog}
  fi
else
  cp ${nfaln} ./
  localn=${nfrad1}
fi

# test presence of SSE3/AVX/AVX2 instruction support
if [[ ! -z "$(grep -o avx2 /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-AVX2)" ]] ; then
  raxmlflav='-AVX2 -U'
elif [[ ! -z "$(grep -o avx /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-AVX)" ]] ; then
  raxmlflav='-AVX -U'
elif [[ ! -z "$(grep -o sse3 /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-SSE3)" ]] ; then
  raxmlflav='-SSE3 -U'
else
  raxmlflav=''
fi

raxmlbin="raxmlHPC-PTHREADS${raxmlflav} -T ${nbthreads}"

raxmloptions="-n ${nfrad2} -m ${model} -p 1753"

idseqgreppat='exactly identical$'
idseqsedpat='s/IMPORTANT WARNING: Sequences \(.\+\) and \(.\+\) are exactly identical/\1\t\2/g'
if [[ "${reducedaln}" == "true" ]] ; then
  ## reduce the alignement and record which sequences were seen as duplicates
  raxmlcall0="${raxmlbin} -s ${localn} ${raxmloptions} -f c && grep '${idseqgreppat}' RAxML_info.${nfrad2} | sed -e '${idseqsedpat}' > RAxML_identical_sequences.${nfrad2}"
else
  raxmlcall0="# NOT reducing alignment to unique sequences before start"
  raxmlbin="${raxmlbin} --silent"
fi
## search for global ML tree
raxmlcall1="${raxmlbin} -s ${localn} ${raxmloptions}"
# search for {rapid|parametric} bootstrap trees
raxmlcall2="${raxmlbin} -s ${localn} ${raxmloptions} -${bootstrapalgo} 987987 -N 100"
## map bootstraps on ML tree
raxmlcall3="${raxmlbin} -s ${localn} ${raxmloptions} -f b -z RAxML_bootstrap.${nfrad2} -t RAxML_bestTree.${nfrad2}"
## root ML tree
raxmlcall4="${raxmlbin} -s ${localn} ${raxmloptions} -f I -t RAxML_bipartitionsBranchLabels.${nfrad2}"
## generic end
raxmlcallz=""

### pipeline
raxmlcalls=(raxmlcall0 raxmlcall1 raxmlcall2 raxmlcall3 raxmlcall4 raxmlcallz)

smallfam=''
status=${?}
for i in {0..5} ; do
  if [ ${status} != 0 ] ; then
    echo "!!! ERROR : during former RAxML call ; exit now"
    exit 1
  else  
    if [ -e "RAxML_info.${nfrad2}" ] ; then 
      mv -f RAxML_info.${nfrad2} RAxML_info.${i}.${nfrad2} 
    fi
    eval raxmlcall='$'${raxmlcalls[$i]}
    echo ""
    echo "#####"
    echo "# RAxML call: ${raxmlcall}"
    eval ${raxmlcall}
    status=${?}
    if [ ${i} -eq 0 ] ; then
      ## change alignment to be further considered by RAxML as the reduced one: !!! the resulting tree will NOT contain any duplicates !!!
      if [[ -e ${localn}.reduced ]] ; then
        nbnrseq=$(head -n1 ${localn}.reduced | cut -d' ' -f1)
        if [[ $nbnrseq -lt 4 ]] ; then
          echo "WARNING: Reduced alignment is too small to pass further steps ; copy the identical_sequences file to '${outputdir}/identical_sequences/' and stops here"
          mv -f ${outputdir}/bulk/RAxML_identical_sequences.${nfrad2} ${outputdir}/identical_sequences/
          smallfam='true'
          break
        else
          echo "# Found $(wc -l RAxML_identical_sequences.${nfrad2} | cut -d':' -f2) redundant sequences; replace input alignment '${localn}' by '${localn}.reduced'."
          repllocaln="mv ${localn} ${localn}.full ; mv ${localn}.reduced ${localn}"
          echo ${repllocaln}
          eval ${repllocaln}
        fi
      fi
    fi
  fi
done &>> ${raxlog}

if [ -z ${smallfam} ] ; then
  for resulttag in ${resulttags[@]} ; do
    mv -f ${outputdir}/bulk/RAxML_${resulttag}.${nfrad2} ${outputdir}/${resulttag}/
  done
else
  echo ${nfrad2} >> ${outputdir}/bulk/reduced_alignment_is_too_small_fams
fi

### end loop over alignments
done
