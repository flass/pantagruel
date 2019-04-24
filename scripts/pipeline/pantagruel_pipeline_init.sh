#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

#### Inherit mandatory environment variables / parameters from parent process = patagruel master script

# derive other important environmnet variables
export ptgscripts="${ptgrepo}/scripts"
source ${ptgscripts}/pipeline/pantagruel_environ_defaults.sh

if [ -z ${initfile} ] ; then 
  templateenv=${ptgscripts}/pipeline/environ_pantagruel_template.sh
else
  templateenv=${initfile}
fi


#~ export PYTHONPATH=$PYTHONPATH:"${ptgrepo}/python_libs"
#~ cd ${ptgrepo} ; export ptgversion=$(git log | grep commit) ; cd - > /dev/null # ptgversion variable should be inherited from master script
# create head folders
export ptgdb=${ptgroot}/${ptgdbname}
export ptglogs=${ptgdb}/logs
export ptgtmp=${ptgdb}/tmp
mkdir -p ${ptgdb}/ ${ptglogs}/ ${ptgtmp}/

echo 'here'
echo $cladesuppdef
echo $collapseCladeParams

# evaluate 'freely specified' collapsing parameters; expect some define paatern though
# e.g.:  eval 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'
for param in 'cladesupp' 'subcladesupp' 'criterion' 'withinfun' ; do
  # parse free text to only evaluate code defining 
  parameqval=$(python << EOF
parameqvals="$collapseCladeParams".lower().split()
for parameqval in parameqvals:
  pareqval = parameqval.strip(';,\t\n ')
  if '=' not in pareqval: continue
  par, val = pareqval.split('=')
  if par=='$param':
    print '%s=%s'%(par, val)
    break
EOF
)
  if [ ! -z "${parameqval}" ] ; then
    echo "export ${parameqval}"
    eval "export ${parameqval}"
  else
    eval "paramdef=\${${param}def}"
    echo "export ${param}=${paramdef}"
    eval "export ${param}=${paramdef}"
  fi
done

#### Set facultative environment variables / parameters
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
echo "## Pantagruel database '${ptgdbname}'" > ${envsourcescript}
echo "## built with Pantagruel version '${ptgversion}'; source code available at 'https://github.com/flass/pantagruel'" >> ${envsourcescript}

ptgversinit="${ptgversion}"

rm -f ${ptgtmp}/sedenvvar.sh
echo -n "cat ${templateenv}" > ${ptgtmp}/sedenvvar.sh
for var in 'ptginitcmd' 'ptgdbname' 'ptgroot' 'ptgrepo' 'ptgversinit' 'myemail' 'famprefix' \
     'ncbitax' 'ncbiass' 'listncbiass' 'customassemb' 'refass' 'listrefass' \
     'pseudocoremingenomes' 'userreftree' 'coreseqtype' 'poplgthresh' 'poplgleafmul' 'popbsthresh' 'rootingmethod' \
     'chaintype' 'cladesupp' 'subcladesupp' 'criterion' 'withinfun' 'hpcremoteptgroot' ; do
echo -n " | sed -e \"s#REPLACE${var}#\${${var}}#\"" >> ${ptgtmp}/sedenvvar.sh
done
echo -n " >> ${envsourcescript}" >> ${ptgtmp}/sedenvvar.sh
bash < ${ptgtmp}/sedenvvar.sh
## load generic environment variables derived from the above
source ${envsourcescript}
if [ -z ${initfile} ] ; then 
  echo "created init file at '${envsourcescript}'"
else
  echo "updated init file save at '${envsourcescript}'"
fi

## check genome input
if [ ! -z "${ncbiass}" ] ; then
  # detected input RefSeq genome folder
    for dass in $(ls -A ${ncbiass}/) ; do
    # check validity of the folder file structure
      if [ ! -d ${ncbiass}/${dass} ] ; then
        echo "Error: RefSeq assemblies must be provided as folders including sequence and annotation files;"
        echo " '${ncbiass}/${dass}' is not a directory; exit now"
        exit 1
      fi
      for ext in genomic.fna genomic.gbff genomic.gff cds_from_genomic protein.faa ; do
       if [ -z "$(ls -A ${ncbiass}/${dass}/${dass}_${ext}* 2> /dev/null)" ] ; then
        echo "Error: could not detect file named like '${dass}_${ext}*' extension in custom assembly folder '${ncbiass}/${dass}/';"
        echo " Input RefSeq assemblies must include a file with extension *_${ext}[.gz] ; exit now"
        exit 1
       fi
      done
    done
fi
if [ ! -z "${listncbiass}" ] ; then
  rm -f ${listncbiass}_wrongFormatAssIds
  grep -v "GC[AF]_[0-9]\{9,\}" ${listncbiass} > ${listncbiass}_wrongFormatAssIds
  if [ -s ${listncbiass}_wrongFormatAssIds ] ; then
    echo "the list of NCBI Assembly accession ids provided through option -L has uncorrectly formatted entries:"
    cat ${listncbiass}_wrongFormatAssIds
    echo "will not be able to download these from NCBI FTP site; exit now"
    exit 1
  fi
fi
if [ ! -z "${customassemb}" ] ; then
  # detected input custom genome folder
  # check validity of the folder file structure
  if [ ! -d "${customassemb}/contigs" ] ; then
    echo "Error: custom assembly folder '${customassemb}/' does not contain mandatory folder 'contigs/'"
    echo " custom assemblies must always provide contigs"
    echo " please copy/link raw sequence (in multi-fasta format) files of custom (user-generated) assemblies into '${customassemb}/contigs/'"
    echo " exit now"
    exit 1
  else
    if [ ! -e "${straininfo}" ] ; then
     echo "Error: '${straininfo}' is missing"
     echo " assembly_id,genus,species,strain,taxid,locus_tag_prefix" | tr ',' '\t' > ${straininfo}
     echo " and fill up the information table ${straininfo} (tab-delimited fields) according to header:"
     cat ${straininfo}
     if [[ "$(ls -A "${customassemb}/contigs/" 2>/dev/null)" ]] ; then
      for allcontigs in `ls ${customassemb}/contigs/` ; do
        gproject=${allcontigs%%.fa*}
        echo "${gproject}\t\t\t\t\t" >> ${straininfo}
      done
      echo " a tab-delimited template was prepared in file '${straininfo}' with assembly_id values copied from names of files found in '${customassemb}/contigs/'"
     fi
     echo " exit now"
     exit 1
    fi
  fi
  if [ -d "${customassemb}/annotations" ] ; then
    for dass in $(ls -A ${customassemb}/annotations/) ; do
      if [ ! -d ${customassemb}/annotations/${dass} ] ; then
        echo "Annotations for custom assemblies must be provided as folders of annotation files;"
        echo " '${customassemb}/annotations/${dass}' is not a directory; exit now"
        exit 1
      fi
      if [ -z "$(ls -A ${customassemb}/annotations/${dass}/*.gff 2> /dev/null)" ] ; then
        echo "Error: could not detect GFF file (with .gff extension) in custom assembly folder '${customassemb}/annotations/${dass}/';"
        echo " Annotations for custom assemblies must include at least a GFF file ; exit now"
        exit 1
      fi
    done
  fi
fi
