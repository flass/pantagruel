#!/usr/bin/env bash

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
source ${ptgscripts}/pipeline/environ_pantagruel_defaults.sh
source ${ptgscripts}/pipeline/pantagruel_pipeline_functions.sh

if [ -z ${initfile} ] ; then 
  templateenv=${ptgscripts}/pipeline/environ_pantagruel_template.sh
else
  echo "Use user-provided template for config  file: ${initfile}"
  templateenv=${initfile}
fi

# create head folders
export ptgdb=${ptgroot}/${ptgdbname}
if [ "${runmode}" != 'wasrefresh' ] ; then
  checkfoldersafe ${ptgdb}
fi


export ptglogs=${ptgdb}/logs
export ptgtmp=${ptgdb}/tmp
mkdir -p ${ptgdb}/ ${ptglogs}/ ${ptgtmp}/

# evaluate 'freely specified' collapsing parameters; expect some define paatern though
# e.g.:  eval 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'
for param in 'cladesupp' 'subcladesupp' 'criterion' 'withinfun' ; do
  # parse free text to only evaluate code defining 
  parameqval=$(python2.7 << EOF
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
    [[ "${runmode}" != 'wasrefresh' && "${chaintype}" == 'collapsed' ]] && echo "${parameqval}"
    eval "export ${parameqval}"
  else
    eval "paramdef=\${${param}def}"
    [[ "${runmode}" != 'wasrefresh' && "${chaintype}" == 'collapsed' ]] && echo "${param}=${paramdef}"
    eval "export ${param}=${paramdef}"
  fi
done

#### Set facultative environment variables / parameters
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh

head -n 1 ${templateenv} > ${envsourcescript}
echo "## Pantagruel database '${ptgdbname}'" >> ${envsourcescript}
echo "## built with Pantagruel version '${ptgversion}'; source code available at 'https://github.com/flass/pantagruel'" >> ${envsourcescript}

export ptgversinit="${ptgversion}"

rm -f ${ptgtmp}/sedenvvar.sh
echo -n "tail -n +16 ${templateenv}" > ${ptgtmp}/sedenvvar.sh
for var in 'ptginitcmd' 'ptgdbname' 'ptgroot' 'ptgrepo' 'ptgversinit' 'myemail' 'famprefix' 'pathtoipscan' \
     'ncbitax' 'ncbiass' 'listncbiass' 'customassemb' 'refass' 'listrefass' 'updatedbfrom' 'customstraininfo' \
     'pseudocoremingenomes' 'userreftree' 'coreseqtype' 'poplgthresh' 'poplgleafmul' 'popbsthresh' 'rootingmethod' 'snpali' \
     'chaintype' 'cladesupp' 'subcladesupp' 'criterion' 'withinfun' 'hpcremoteptgroot' 'genefamlist' 'preferredgenomes' 'recmethod' 'maxreftreeheight' ; do
echo -n " | sed -e \"s#REPLACE${var}#\${${var}}#\"" >> ${ptgtmp}/sedenvvar.sh
done
echo -n " >> ${envsourcescript}" >> ${ptgtmp}/sedenvvar.sh
bash < ${ptgtmp}/sedenvvar.sh

if [ ! -z ${extravars} ] ; then
  # add custom env variables
  echo "custom environment variables:"
  echo "${extravars}" | tr ',' '\n' | while read exp ; do
    if [[ -z "$(echo ${exp} | grep '=')" || ! -z "$(echo ${exp} | grep '^export ')" ]] ; then 
	  echo "Error: '${exp}' is not a properly formatted variable assignment expression; also please discard any 'export' statement at the begining. Exit now."
	  exit 1
	fi
    echo "export ${exp}"
    echo "export ${exp}" >> ${envsourcescript}
  done
fi

## load generic environment variables derived from the above
source ${envsourcescript}
checkexec "Succesfully created config file at '${envsourcescript}'"

if [ "${runmode}" == 'wasrefresh' ] ; then
  echo "refresh mode: skip input folder checks."
  exit 0
elif [[ "${runmode}" == 'force' ]] ; then
  echo "FORCE mode: pantagruel init will skip the checks on input dataset integrity; the following tasks may fail if these tests are not clear."
  exit 0
fi

## check genome input
for srcass in "${ncbiass}" "${refass}" ; do
if [ ! -z "${srcass}" ] ; then
  inass="$(ls -A "${srcass}/" 2>/dev/null)"
  if [ -z "${inass}" ] ; then
    echo "Error: the folder '${srcass}' provided for custom reference genome assembly set (through option '--refseq_ass4annot') is empty; pease provide "
  fi
  # detected an input RefSeq genome folder
    for dass in ${inass} ; do
      # check validity of the folder file structure
      if [[ "${dass}" == 'all_assemblies_organism_names' || "${dass}" == 'genome_assemblies_list' ]] ; then
        echo "found file '${dass}' presumably generated from a previous Pantagruel run; this file will be deleted now and regenrated later"
        rm ${srcass}/${dass}
        continue
      elif [ ! -d ${srcass}/${dass} ] ; then
        echo "Error: RefSeq assemblies must be provided as folders including sequence and annotation files;"
        echo " '${srcass}/${dass}' is not a directory; exit now"
        exit 1
      fi
      for ext in genomic.fna genomic.gbff genomic.gff cds_from_genomic protein.faa ; do
       if [ -z "$(ls -A ${srcass}/${dass}/${dass}_${ext}* 2> /dev/null)" ] ; then
        echo "Error: could not detect file named like '${dass}_${ext}*' extension in assembly folder '${srcass}/${dass}/';"
        echo " Input RefSeq assemblies must include a file with extension *_${ext}[.gz] ; exit now"
        exit 1
       fi
      done
    done
fi
done

for listsrcass in "${listncbiass}" "${listrefass}" ; do
if [ ! -z "${listsrcass}" ] ; then
  wfalist=${ptgtmp}/$(basename ${listsrcass})_wrongFormatAssIds
  rm -f ${wfalist}
  grep -v "GC[AF]_[0-9]\{9,\}" ${listsrcass} > ${wfalist}
  if [ -s ${wfalist} ] ; then
    echo "the list of NCBI Assembly accession ids provided through option -L has uncorrectly formatted entries:"
    cat ${wfalist}
    echo "will not be able to download these from NCBI FTP site; exit now"
    exit 1
  fi
fi
done

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
     echo -e "assembly_id\tgenus\tspecies\tstrain\ttaxid\tlocus_tag_prefix" > ${straininfo}
     echo " a template information table '${straininfo}' was created"
     if [[ "$(ls -A "${customassemb}/contigs/" 2>/dev/null)" ]] ; then
      for allcontigs in `ls ${customassemb}/contigs/` ; do
        gproject=${allcontigs%%.fa*}
        echo -e "${gproject}\t\t\t\t\t" >> ${straininfo}
      done
     fi
     echo " with assembly_id values copied from names of files found in '${customassemb}/contigs/';"
     echo " please fill it up according to header (and respecting already present tabs delimiting fields):"
     cat ${straininfo}
     echo " exit now"
     exit 1
    else
      # check that proposed locus-tag prefix are free from characters '-' and '_'
      python2.7 << EOF
nfin = '${straininfo}'
errprefix = "Error in format of strain info file '%s': "%nfin
errsuffix = "\nPlease edit the file accordingly."
dstrains = {}
expectedfields = set(['assembly_id','genus','species','strain','taxid','locus_tag_prefix'])
with open(nfin, 'r') as fin:
  header = fin.readline().rstrip('\n').split('\t')
  missingfields = expectedfields - set(header)
  extrafields = set(header) - expectedfields
  if missingfields:
    sextra = "\nExtra fields were detected: %s. Maybe the header was misspelled?"%repr(list(extrafields)) if extrafields else ''
    raise ValueError, "%s\nthe fields %s are missing from the header.%s%s"%(errprefix, repr(list(missingfields)), sextra, errsuffix)
  ltpi = header.index('locus_tag_prefix')
  for line in fin:
    lsp = line.rstrip('\n').split('\t')
    ltp = lsp[ltpi]
    if ('-' in ltp) or ('_' in ltp):
      raise ValueError, "%s the characters '-' and '_' are forbiden in the 'locus_tag_prefix' field. %s"%(errprefix, errsuffix)

EOF
    checkexec "Custom strain info file detected; Error in format" "Custom strain info file detected; format validated"
    fi
  fi
  if [ -d "${custannot}" ] ; then
    for dass in $(ls -A ${custannot}/) ; do
      if [ ! -d ${custannot}/${dass} ] ; then
        echo "Annotations for custom assemblies must be provided as folders of annotation files;"
        echo " '${custannot}/${dass}' is not a directory; exit now"
        exit 1
      fi
      if [ -z "$(ls -A ${custannot}/${dass}/*.gff 2> /dev/null)" ] ; then
        echo "Error: could not detect GFF file (with .gff extension) in custom assembly folder '${custannot}/${dass}/';"
        echo " Annotations for custom assemblies must include at least a GFF file ; exit now"
        exit 1
      fi
    done
  fi
fi
