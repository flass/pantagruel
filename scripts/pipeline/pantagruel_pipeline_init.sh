#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

#### Inherit mandatory environment variables / parameters from parent process = patagruel master script
#~ #### Set mandatory environment variables / parameters
#~ export ptgdbname="${1}"  # database anme (will notably be the name of the top folder)
#~ export ptgroot="${2}"    # source folder where to create the database
#~ export ptgrepo="${3}"    # path to the pantagruel git repository
#~ export myemail="${4}"    # user identity
#~ export famprefix="${5}"  # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.          
#~ export ncbiass="${6}"
#~ export ncbitax="${7}"
#~ export customassemb="${8}"
#~ export refass="${9}"
#~ export chaintype="${10}"
#~ export pseudocoremingenomes="${11}"
#~ export hpcremoteptgroot="${12}"
#~ export collapseCladeParams="${13}"
#~ export userreftree="${14}"
#~ export coreseqtype="${15}"
#~ export poplgthresh="${16}"
#~ export popbsthresh="${17}"
#~ export initfile="${18}"

# derive other important environmnet variables
export ptgscripts="${ptgrepo}/scripts"

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

#### Set facultative environment variables / parameters
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
echo "## Pantagruel database '${ptgdbname}'" > ${envsourcescript}
echo "## built with Pantagruel version '${ptgversion}'; source code available at 'https://github.com/flass/pantagruel'" >> ${envsourcescript}

rm -f ${ptgtmp}/sedenvvar.sh
echo -n "cat ${templateenv}" > ${ptgtmp}/sedenvvar.sh
for var in ptgdbname ptgroot ptgrepo myemail famprefix \
     ncbiass ncbitax customassemb refass chaintype \
     pseudocoremingenomes hpcremoteptgroot collapseCladeParams userreftree coreseqtype \
     poplgthresh popbsthresh ; do
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

# folders for optional custom genomes
if [ -d "${customassemb}/contigs" ] ; then
  if [ ! -e ${straininfo} ] ; then
   echo "assembly_id,genus,species,strain,taxid,locus_tag_prefix" | tr ',' '\t' > ${straininfo}
   echo "please copy/link raw sequence (in multi-fasta format) files of custom (user-generated) assemblies into '${customassemb}/contigs/'"
   echo "and fill up the information table ${straininfo} (tab-delimited fields) according to header:"
   cat ${straininfo}
   if [[ "$(ls -A "${customassemb}/contigs/" 2>/dev/null)" ]] ; then
    for allcontigs in `ls ${customassemb}/contigs/` ; do
      gproject=${allcontigs%%.fa*}
      echo "${gproject}\t\t\t\t\t" >> ${straininfo}
    done
    echo "prepared tab-delimited rows in file '${straininfo}' from files found in '${customassemb}/contigs/'"
   fi
  fi
else
  if [ ! -d "${customassemb}/annotations" ] ; then
   echo "custom assembly folder '${customassemb}/' contains the facultative folder 'annotations/' but not mandatory fodler 'contigs/'"
   echo "Error: annotations mus always be accompanied of contigs ; exit now"
   exit 1
  fi
fi
