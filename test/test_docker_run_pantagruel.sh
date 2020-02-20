#!/bin/bash
cwd=${1}
ptgroot=${2}
ptgdbname=${3}
dockerimage=${4}
ptgrepo=${5}
[ -z ${cwd} ] && cwd=${PWD}									# working directory; defaults to current folder
[ -z ${ptgroot} ] && ptgroot="${cwd}/pantagruel_databases"	# target folder for new database; defaults to ./pantagruel_databases
[ -z ${ptgdbname} ] && ptgdbname="testPTGdatabase"			# name for new database;
#															  defaults to 'testPTGdatabase', which will trigger test of software
[ -z ${dockerimage} ] && dockerimage="panta"				# name of the docker image to use; defaults to 'panta'
[ -z ${ptgrepo} ] && ptgrepo="/pantagruel"					# pantagruel git repository folder containing the pantagruel executables and scripts; 
#															  defaults to the docker image's built-in repo in /pantagruel

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

dockerrunptg="docker run --user ${UID}:${UID} -v ${cwd}:${cwd} -w ${cwd} ${dockerimage}"

${dockerrunptg} ${ptgexe} -r ${ptgroot} -d ${ptgdbname} ${initargs} init
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} --refresh ${refreshargs} init
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 00
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 01
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 02
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 03
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 05
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 06
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 07
${dockerrunptg} ${ptgexe} -i ${ptgdbconfig} 08
