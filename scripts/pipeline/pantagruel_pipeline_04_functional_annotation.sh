#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 26 November 2018

if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source $envsourcescript

cd ${ptgrepo} ; export ptgversion=$(git log | grep commit | cut -d' ' -f2) ; cd -

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
lastIPversion=$(lftp -c "open ${iphost} ; ls -tr ${iploc} ; quit" | tail -n 1 | awk '{print $NF}')
if [ "${currIPversion}" != "${lastIPversion}" ] ; then
  echo "Error: the installed verison of InterProScan (found at $(ls -l `which interproscan` | awk '{print $NF}')) is ${currIPversion}, different from the version currently operated by the EBI: ${lastIPversion}."
  echo "Using this outdated version would cause the look-up service not to work and thus a significant loss of efficiency."
  echo "Please install the most recent version by running again the script '${ptgrepo}/install_dependencies.sh'."
  echo "Exiting now."
  exit 1
fi

export interpro=${funcannot}/InterProScan_${IPversion}
mkdir -p ${interpro}/ ${enttmp}/interpro/ ${entlogs}/interpro/


# segment the entire proteome in batches of reasonable size (as all results are kept in memory before bing written)
mkdir -p ${interpro}/all_complete_proteomes
ln -s ${nrfaacomplete} ${interpro}/
python << EOF
import os
nfinfaa = "${nrfaacomplete}"
outd = "${interpro}/all_complete_proteomes"
nfoutrad = os.path.basename(nfinfaa)
maxnprot = 10000
k = 0
n = 0
finfaa = open(nfinfaa, 'r')
fout = open(os.path.join(outd, nfoutrad)+'.%d'%k, 'w')
for line in finfaa:
  if line.startswith('>'):
      if n > maxnprot:
        fout.close()
        k += 1
        fout = open(os.path.join(outd, nfoutrad)+'.%d'%k, 'w')
        n = 0
      n += 1
  fout.write(line)

fout.close()
EOF

for nrfaa in `ls ${interpro}/all_complete_proteomes/*.faa*` ; do
 interproscan -T ${enttmp}/interpro -i ${nrfaa} -b ${nrfaa} --formats TSV --iprlookup --goterms --pathways &> ${entlogs}/interpro/$(basename ${nrfaa})_interpro.log
done &


## load the annotation data into the database
#~ python ${ptgscripts}/pantagruel_sqlitedb_load_interproscan_annotations.py ${sqldb} ${interpro} ${IPversion}
python ${ptgscripts}/pantagruel_load_interproscan_annotations.py --sqlite_db ${sqldb} --ipscan_annot_files "${interpro}/all_complete_proteomes/*.tsv" --ipscan_version ${IPversion}

