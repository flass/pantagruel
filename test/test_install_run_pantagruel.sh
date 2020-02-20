#!/bin/bash

cwd=${1}
ptgroot=${2}
ptgdbname=${3}
ptgrepodest=${4}
[ -z ${cwd} ] && cwd=${PWD}									# working directory; defaults to current folder
[ -z ${ptgroot} ] && ptgroot="${cwd}/pantagruel_databases"	# target folder for new database; defaults to ./pantagruel_databases
[ -z ${ptgdbname} ] && ptgdbname="testPTGdatabase"			# name for new database;
#															  defaults to 'testPTGdatabase', which will trigger test of software
[ -z ${ptgrepodest} ] && ptgrepodest=${PWD}	# where to clone the pantagruel git repository folder containing the pantagruel executables and scripts; 
#                                             defaults to the current folder
ptgrepo="${ptgrepodest}/pantagruel"

ptgdbconfig="${ptgroot}/${ptgdbname}/environ_pantagruel_${ptgdbname}.sh"

if [ ${ptgdbname} == "testPTGdatabase" ] ; then
  echo "will run the test of software within the Docker container"
  # location of Pantagruel executable and test data
  ptgdata="${ptgrepo}/data"
  initargs="-L ${ptgdata}/NCBI_Assembly_accession_ids_test_10Brady -a ${ptgdata}/custom_genomes"
  refreshargs="-c -g ${ptgdata}/test_genefam_list -q 0.75"
else
  # use this to run on actual data
  initargs="${@}"  # all trailing script arguments
fi

echo "test installing pipeline"
mkdir -p ${ptgrepodest}
cd ${ptgrepodest}
git clone https://github.com/flass/pantagruel.git
cd -
${ptgrepo}/install_dependencies.sh pantagruel_pipeline
echo "test whether the pantagruel command is available from the PATH"
pantagruel -h
if [ ${?} -gt 0 ] ; then
  echo "WARNING: pantagruel command is NOT AVAIALBE from the PATH"
  ptgexe="${ptgrepo}/pantagruel"
else
  echo "pantagruel command is available from the PATH"
  ptgexe="pantagruel"
fi

echo "test running pipeline"
${ptgexe} -d testPTGdatabase -r ./ -f PANTAGFAM -I e.mail@institu.ti.on -L ./pantagruel_pipeline/pantagruel/data/NCBI_Assembly_accession_ids_test_10Brady -a ./pantagruel_pipeline/pantagruel/data/custom_genomes init
${ptgexe} -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh all
