#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 26 November 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

checkptgversion
checkfoldersafe ${funcannot}

if [ -z "${interpro}" ] ; then
  echo "InterproScan was not found on this machine: cannot run this (facultative) task of Pantagruel pipeline; exit now"
  exit 1
fi

# make sure to use the correct Java version (incompatibility should have been detected during install of Pantagruel)
if [ ! -z ${JAVA4INTERPROSCAN} ] ; then
  # if default version of java is not the one expected by InterProScan, a line should have been added to ~/.bashrc to define the JAVA4INTERPROSCAN environment variable
  export PATH=$(dirname ${JAVA4INTERPROSCAN}):$PATH
fi

############################
## 04. Functional annotation
############################

mkdir -p ${funcannot}/
export currIPversion=$(interproscan --version | head -n 1 | sed -e 's/InterProScan version //')
if [ -z "$currIPversion" ] ; then
  echo "Error: unable to dertermine version of Interproscan; please verify the program is correctly installed ; exiting now."
  exit 1
fi
iphost="ftp://ftp.ebi.ac.uk"
iploc="pub/software/unix/iprscan/5/"
lastIPversion=$(lftp -c "open ${iphost} ; ls -tr ${iploc}[0-9].[0-9][0-9]*-[0-9][0-9]*.[0-9] ; quit" | tail -n 1 | awk '{print $NF}')
if [ "${currIPversion}" != "${lastIPversion}" ] ; then
  echo "Error: the installed verison of InterProScan (found at $(ls -l `which interproscan` | awk '{print $NF}')) is '${currIPversion}', different from the version currently published on the EBI FTP: '${lastIPversion}'."
  echo "Using this outdated version would cause the look-up service not to work and thus a significant loss of efficiency."
  echo "Please install the most recent version by running again the script '${ptgrepo}/install_dependencies.sh'."
  echo "Exiting now."
  exit 1
fi

export interpro=${funcannot}/InterProScan_${IPversion}
mkdir -p ${interpro}/ ${ptgtmp}/interpro/ ${ptglogs}/interpro/


# segment the entire proteome in batches of reasonable size (as all results are kept in memory before being written)
mkdir -p ${interpro}/all_complete_proteomes/
ln -s $(realpath --relative-to=${interpro} ${allfaarad}.nr.faa) ${interpro}/
python << EOF
import os
nfinfaa = "${allfaarad}.nr.faa"
outd = "${interpro}/all_complete_proteomes"
nfoutrad = os.path.basename(nfinfaa)
maxnprot = 10000
k = 0
n = 0
finfaa = open(nfinfaa, 'r')
fout = open(os.path.join(outd, nfoutrad)+'.%d'%k, 'w')
for line in finfaa:
  if line.startswith('>'):
      if n >= maxnprot:
        fout.close()
        k += 1
        fout = open(os.path.join(outd, nfoutrad)+'.%d'%k, 'w')
        n = 0
      n += 1
  fout.write(line)

fout.close()
EOF

if [ ! -z "${ptgthreads}" ] ; then
  ipcpu="--cpu ${ptgthreads}"
fi
for nrfaa in $(ls ${interpro}/all_complete_proteomes/*.faa* | grep -v 'tsv') ; do
 if [[ "${resumetask}" == 'true' && -s "${nrfaa}.tsv" ]] ; then
   echo "Found the output of a previous InterProScan run for protein flat file ${nrfaa}: ${nrfaa}.tsv ; skip computation."
 else
   interproscan -T ${ptgtmp}/interpro -i ${nrfaa} -b ${nrfaa} --formats TSV ${ipcpu} --iprlookup --goterms --pathways &> ${ptglogs}/interpro/$(basename ${nrfaa})_interpro.log
 fi
done


## load the annotation data into the database
#~ python ${ptgscripts}/pantagruel_sqlitedb_load_interproscan_annotations.py ${sqldb} ${interpro} ${IPversion}
python ${ptgscripts}/pantagruel_load_interproscan_annotations.py --sqlite_db ${sqldb} --ipscan_annot_files "${interpro}/all_complete_proteomes/*.tsv" --ipscan_version ${IPversion}

