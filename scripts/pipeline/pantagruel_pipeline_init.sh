#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

#### Set mandatory environment variables / parameters
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgrepo="$3"    # path to the pantagruel git repository
export myemail="$4"    # user identity
export famprefix="$5"  # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.          
export ncbiass="$6"
export ncbitax="$7"
export customassemb="$8"
export initfile="$9"

# derive other important environmnet variables
export ptgscripts="${ptgrepo}/scripts"

if [ -z ${initfile} ] ; then 
  templateenv=${ptgscripts}/pipeline/environ_pantagruel_template.sh
else
  templateenv=${initfile}
fi


#~ export PYTHONPATH=$PYTHONPATH:"${ptgrepo}/python_libs"
cd ${ptgrepo} ; export ptgversion=$(git log | grep commit) ; cd - > /dev/null
# create head folders
export ptgdb=${ptgroot}/${ptgdbname}
export ptglogs=${ptgdb}/logs
export ptgtmp=${ptgdb}/tmp
mkdir -p ${ptgdb}/ ${ptglogs}/ ${ptgtmp}/

#### Set facultative environment variables / parameters
export pseudocoremingenomes=''       # defaults to empty variable in which case will be set INTERACTIVELY at stage 04.core_genome of the pipeline
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh

rm -f ${ptgtmp}/sedenvvar.sh
echo -n "cat ${templateenv}" > ${ptgtmp}/sedenvvar.sh
for var in ptgdbname ptgroot ptgrepo myemail famprefix customassemb ncbiass ncbitax ; do
echo -n " | sed -e \"s#REPLACE${var}#\${${var}}#\"" >> ${ptgtmp}/sedenvvar.sh
done
echo -n " > ${envsourcescript}" >> ${ptgtmp}/sedenvvar.sh
bash < ${ptgtmp}/sedenvvar.sh
## load generic environment variables derived from the above
source ${envsourcescript}
if [ -z ${initfile} ] ; then 
  echo "created init file at '${envsourcescript}'"
else
  echo "updated init file save at '${envsourcescript}'"
fi

# folders for optional custom genomes
mkdir -p ${customassemb}/contigs/
if [ ! -e ${straininfo} ] ; then
 echo "sequencing.project.id,genus,species,strain,taxid,locus_tag_prefix" | tr ',' '\t' > ${straininfo}
 echo "please copy/link raw sequence (in multi-fasta format) files of custom (user-generated) assemblies into '${customassemb}/contigs/'"
 echo "and fill up the information table ${straininfo} (tab-delimited fields) according to header:"
 cat ${straininfo}
fi

