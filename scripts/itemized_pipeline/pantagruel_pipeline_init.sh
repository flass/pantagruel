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
export myemail="me.myself@respectable-institu.ti.on"
export raproot="/path/to/base/folder/for/all/this/business"
export ptgrepo="/path/to/pantagruel_repository"
export famprefix="REPLACEfamprefix"
export rapdbname="aBoldDatabaseName" # mostly name of the top folder
export famprefix="ABCDE"             # alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs.

# derive other important environmnet variables
export ptgscripts="${ptgrepo}/scripts"
export PYTHONPATH=$PYTHONPATH:${ptgscripts}
cd ${ptgrepo} ; export ptgversion=$(git log | grep commit) ; cd -

# create head folders
export rapdb=${raproot}/${rapdbname}
export raplogs=${rapdb}/logs
export raptmp=${rapdb}/tmp
mkdir -p ${rapdb}/ ${raplogs}/ ${raptmp}/

#### Set facultative environment variables / parameters
export pseudocoremingenomes=''       # defaults to empty variable in which case will be set INTERACTIVELY at stage 04.core_genome of the pipeline
envsourcescript=${raproot}/environ_pantagruel.sh

rm -f ${raptmp}/sedenvvar.sh
echo -n "cat ${ptgscripts}/environ_pantagruel_template.sh" > ${raptmp}/sedenvvar.sh
for var in myemail raproot ptgscripts famprefix rapdbname famprefix pseudocoremingenomes ; do
echo -n " | sed -e \"s#REPLACE${var}#${var}#\"" >> ${raptmp}/sedenvvar.sh
done
echo -n " > ${envsourcescript}" >> ${raptmp}/sedenvvar.sh
bash < ${raptmp}/sedenvvar.sh
## load generic environment variables derived from the above
source ${envsourcescript}
cat "source ${envsourcescript}" >> ~/.bashrc

# folders for optional custom genomes
export customassemb=${raproot}/user_genomes
mkdir -p ${customassemb}/contigs/
echo "sequencing.project.id,genus,species,strain,taxid,locus_tag_prefix" | tr ',' '\t' > ${straininfo}

echo "please copy/link raw sequence (in multi-fasta format) files of custom (user-generated) assemblies into ${customassemb}/contigs/"
echo "and fill up the information table ${straininfo} (tab-delimited fields) according to header:"
cat ${straininfo}
