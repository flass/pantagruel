#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 26 November 2018

if [ -z "$2" ] ; echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source $envsourcescript

cd ${ptgrepo} ; export ptgversion=$(git log | grep commit | cut -d' ' -f2) ; cd -

############################
## 04. Functional annotation
############################

export funcannot=${ptgdb}/04.functional
mkdir -p ${funcannot}/
export IPversion=$(interproscan --version | head -n 1 | sed -e 's/InterProScan version //')
if [ -z "$IPversion" ] ; then
  echo "Error unable to dertermine version of Interproscan; please verify the program is correctly installed ; exiting now."
  exit(1)
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

