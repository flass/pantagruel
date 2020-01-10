#!/bin/bash
cwd=${1}	
ptgrepo=${2}							  
ptgroot=${3}							  
ptgdbname=${3}							  
[ -z ${cwd} ] && cwd=${PWD}									# working directory, defaults to current folder
[ -z ${ptgrepo} ] && ptgrepo="${cwd}/software/pantagruel"	# pantagruel git repository folder, defaults to ./software/pantagruel
[ -z ${ptgroot} ] && ptgroot="${cwd}/pantagruel_databases"	# target folder for new database, defaults to ./pantagruel_databases
[ -z ${ptgdbname} ] && ptgdbname="testPTGdatabase"			# name for new database, defaults to 'testPTGdatabase', which will trigger test of software

ptgexe="${ptgrepo}/pantagruel"
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

dockerrunptg="docker run --user $USER -v $cwd:$cwd -w $cwd pantagruel"

${dockerrunptg} ${ptgexe} -r ${ptgroot} -d ${ptgdbname} ${initargs} init
cd ${ptgrepo}/ && git pull && git submodule update && cd -
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} --refresh ${refreshargs} init
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 00
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 01
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 02
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 03
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 05
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 06
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 07
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 08
