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
## 09. Functional annotation
############################

export funcannot=${ptgdb}/09.functional
mkdir -p ${funcannot}
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
python ${ptgscripts}/pantagruel_sqlitedb_load_interproscan_annotations.py ${sqldb} ${interpro} ${IPversion}
#~ python ${ptgscripts}/pantagruel_load_interproscan_annotations.py --sqlite_db ${sqldb} --ipscan_annot_files "${interpro}/all_complete_proteomes/*.tsv" --ipscan_version ${IPversion}


sqlite3 -cmd ".mode tab" ${sqldb} "select distinct locus_tag, interpro_id, interpro_description, go_terms, pathways 
from coding_sequences 
inner join functional_annotations using (nr_protein_id) 
inner join interpro_terms using (interpro_id) 
where interpro_description like '%arsen%' order by locus_tag, interpro_id;" > ${funcannot}/arsen_interpro_matches_go_terms_pathways.tab

## generate background term distribution for clades i.e. list all genes in their pangenome and associated terms
export goterms=${funcannot}/GeneOntology
mkdir -p ${goterms}
# for the whole dataset
sqlite3 -cmd ".mode tab" ${sqldb} "select distinct locus_tag, go_id from coding_sequences 
left join functional_annotations using (nr_protein_id) 
left join interpro2GO using (interpro_id) ;" > ${goterms}/${ngenomes}-genomes_pangenome_terms.tab
# for clades of the $focus clade 41NeoPseudo
export claderefgodir=${goterms}/${focus}_clade_go_term_reference_sets
mkdir -p ${claderefgodir}/
tail -n +2 ${corefocus}/${focus}_clade_defs | while read cla cladedef siscladedef ; do
  claspeset="'$(echo $cladedef | sed -e "s/,/','/g")'"
  echo $cla $claspeset
  cladest=${claderefgodir}/${focus}_${cla}_pangenome_terms.tab
  q="select distinct locus_tag, go_id from coding_sequences 
  inner join replicons using (genomic_accession)
  inner join assemblies using (assembly_id)
  left join functional_annotations using (nr_protein_id) 
  left join interpro2GO using (interpro_id) 
  where code in (${claspeset})"
  sqlite3 -cmd ".mode tab" ${sqldb} "${q};" > ${cladest}
  sqlite3 -cmd ".mode tab" ${sqldb} "${q} and go_id not null;" > ${cladest}_nonull
  ls -lh ${cladest}
done

## test GO term enrichment in gene sets
catrefile (){
  for cla in $1 ; do 
    ls ${claderefgodir}/${focus}_${cla}_pangenome_terms.tab_nonull
  done | xargs | tr ' ' ','
}
# first test general enrichment of a clade pangenome against the larger pangenome
${ptgscripts}/clade_specific_genes_GOterm_enrichment_test.r \
--study_annots "$( catrefile 'clade8 clade9 clade10')" \
--population_annots "$( catrefile 'clade4')" \
--out ${goterms}/${focus}_clade_go_term_enriched_vs_clade4_pangenome \
--algo "weight01" --stat "Fisher" &> ${raplogs}/clade_specific_genes_GOterm_enrichment_test.log &
#~ [1] 21044

# then refresh clade-specific gene tables, adding GO term and pathway annotations
${ptgscripts}/get_clade_specific_genes.r ${orthomatrad}_genome_counts.no-singletons.mat ${sqldb} ${orthocolid} ${coregenome}/${focus}/${focus} ${orthomatrad}

