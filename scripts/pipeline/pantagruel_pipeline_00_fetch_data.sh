#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"
export ptgroot="$2"
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

#################################
## 00. Data download and grooming
#################################

if [ ! -e ${ncbitax}/names.dmp ] ; then
echo "did not find the relevant taxonomy flat files in '${ncbitax}/'; download the from NCBI Taxonomy FTP"
  ### Download data
  ## using FTP
  # downloading relevant tables from NCBI Taxonomy database 
  user='anonymous'
  pswd=${myemail}
  host='ftp.ncbi.nlm.nih.gov'
  openparam=" -u ${user},${pswd} ${host}"
  source='/pub/taxonomy'
  files='taxcat.tar.gz* taxcat_readme.txt taxdump.tar.gz* taxdump_readme.txt'
  mkdir -p ${ncbitax}
  lftp -c "cd ${source} ; mget -O ${ncbitax}/ ${files}" ${openparam}
  cd ${ncbitax}/
  for tgz in `ls *.tar.gz` ; do md5sum -c ${tgz}.md5 && tar -xzf ${tgz} ; done
  cd -
fi

## using NCBI web interface:
# Search Assembly database using a query defined as convenient:  
# e.g.: 'txid543[Organism:exp] AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter]) AND ("complete genome"[filter])'
# save assemblies (Download Assemblies > Source Database: RefSeq; File type: All file type (including assembly-structure directory))
# as ${ncbiass}/genome_assemblies.tar or as ${ncbiass}/genome_assemblies_DATASETTAG.tar if several datasets have to be merged
# extract assembly data
mkdir -p  ${assemblies}/
cd ${ncbiass}/
# look for assembly folders
ls -A ${ncbiass}/ | grep -v 'assemblies' > ${ncbiass}/genome_assemblies_list
# look for archives containing many assembly folders (datasets downloaded from NCBI Assembly website)
assarch=$(ls genome_assemblies*.tar 2> /dev/null)
if [ "$assarch" ] ; then
  echo "detected archive(s) of genome assemblies (presumably downloaded from NCBI Assembly website):"
  echo "$assarch"
  # extract assembly folders from the archives
  for tarf in `ls genome_assemblies*.tar` ; do 
    tartag=${tarf#genome_assemblies*}
    datasettag=${tartag%.*}
    tar -tf genome_assemblies${datasettag}.tar | cut -d'/' -f2 | sort -u | grep -v "README\|report" > ${ncbiass}/genome_assemblies${datasettag}_list
    tar -xf genome_assemblies${datasettag}.tar 
    mv report.txt report${datasettag}.txt
    assd=(`ls ncbi-genomes-* -d`)
    for dass in `ls ${assd[0]}/ | grep -v README` ; do
     if [ ! -z $dass ] ; then
    rm -rf ${ncbiass}/$dass
    mv ${assd[0]}/$dass ${ncbiass}/
     fi
    done
    rm -r ${assd[0]}/ 
  done
fi

grep "# Organism name" ${ncbiass}/*/*_assembly_stats.txt > ${ncbiass}/all_assemblies_organism_names

# store in centralised folder for NCBI assemblies and just record links
mkdir -p ${assemblies}/
for ass in `cat ${ncbiass}/genome_assemblies*_list` ; do
  ln -s ${ncbiass}/${ass} ${assemblies}/
done


## if the folder of custom/user-provided set of genomes is not empty
if [[ "$(ls -A "${contigs}/" 2>/dev/null)" ]] ; then

  ### uniformly annotate the raw sequence data with Prokka
  # first add all the (representative) proteins in the dataset to the custom reference prot database for Prokka to search for similarities
  ls ${assemblies}/*/*_genomic.gbff.gz > ${indata}/assemblies_genomic_gbffgz_list
  parallel -a ${indata}/assemblies_genomic_gbffgz_list 'gunzip -k'
  # extract protein sequences
  prokka-genbank_to_fasta_db ${assemblies}/*/*_genomic.gbff > ${ptgtmp}/${refgenus}.faa 2> ${ptglogs}/prokka-genbank_to_fasta_db.log
  # cluster similar sequences
  cdhit -i ${ptgtmp}/${refgenus}.faa -o ${ptgtmp}/${refgenus}_representative.faa -T 0 -M 0 -G 1 -s 0.8 -c 0.9 &> ${ptglogs}/cdhit.log
  rm -fv ${ptgtmp}/${refgenus}.faa ${ptgtmp}/${refgenus}_representative.faa.clstr
  # replace database name for genus detection by Prokka 
  prokkabin=$(which prokka)
  prokkablastdb=$(dirname $(dirname $(readlink -f $prokkabin)))/db/genus
  cp -p ${ptgtmp}/${refgenus}_representative.faa ${prokkablastdb}/${refgenus}
  cd ${prokkablastdb}/
  makeblastdb -dbtype prot -in ${refgenus}
  cd -

  # run Prokka 
  mkdir -p ${annot}
  for allcontigs in `ls ${contigs}/` ; do
    gproject=${allcontigs%%.fa*}
    if [ -d ${annot}/${gproject} ] ; then
      echo "found annotation folder '${annot}/${gproject}' ; skip annotation of contigs in '${contigs}/${allcontigs}'"
    else
      echo "will annotate contigs in '${contigs}/${allcontigs}'"
      echo "# $(dateprompt)"
      echo "### assembly: $gproject; contig files from: ${contigs}/${allcontigs}"
      echo "running Prokka..."
      ${ptgscripts}/run_prokka.sh ${gproject} ${contigs}/${allcontigs} ${straininfo} ${annot}/${gproject} ${seqcentre} &> ${ptglogs}/${ptgdbname}_customassembly_annot_prokka.${gproject}.log
      echo "done."
      echo "# $(dateprompt)"
    fi
    if [[ $(grep ${gproject} ${straininfo} | cut -f2- | grep -P '[^\t]' | wc -l) -gt 0 ]] ; then
      echo "fix annotation to integrate region information into GFF files"
      annotgff=($(ls ${annot}/${gproject}/*.gff))
      if [ -z ${annotgff[0]} ] ; then
        echo "Error: missing mandatory GFF file in ${annot}/${gproject}/ folder"
        exit 1
      fi
      python ${ptgscripts}/add_region_feature2prokkaGFF.py ${annotgff[0]} ${annotgff[0]%*.gff}.ptg.gff ${straininfo} ${contigs}/${allcontigs} ${assembler}
      echo "fix annotation to integrate taxid information into GBK files"
      annotfna=($(ls ${annot}/${gproject}/*.fna))
      annotgbk=($(ls ${annot}/${gproject}/*.gbk))
      annotfaa=($(ls ${annot}/${gproject}/*.faa))
      annotffn=($(ls ${annot}/${gproject}/*.ffn))
      if [[ -z "${annotfna[0]}" || -z "${annotgbk[0]}" || -z "${annotffn[0]}" || -z "${annotfaa[0]}" ]] ; then
        echo "At least one of these files is missing in ${annot}/${gproject}/ folder: contig fasta file (.fna), GenBank flat file (gbk), CDS Fasta (ffn) or protein Fasta (faa)."
        echo "Will (re)generate them from the GFF anotation and genomic Fasta sequence; files already present are kept with an added prefix '.original'"
        for annotf in ${annotfna[@]} ${annotgbk[@]} ${annotffn[@]} ${annotfaa[@]}; do
          if [[ ! -z "${annotf}" ]] ; then
            mv ${annotf} ${annotf}.original
          fi
        done
        if [ -z $(grep '##sequence-region' ${annotgff[0]}) ] ; then
		  mv ${annotgff[0]} ${annotgff[0]}.original
          # make GFF source file look more like the output of Prokka
          head -n1 ${annotgff[0]}.original > ${annotgff[0]}
          # add region annotation features
          python << EOF
          fcontig = open('${contigs}/${allcontigs}', 'r')
          outgff = open('${annotgff[0]}', 'a')
          inseq = False
          seqlen = 0
          seqid = None
          for line in fcontig:
            if line.startswith('>'):
              if inseq:
                # write out previous
                outgff.write("##sequence-region %s 1 %d\n"%(seqid, seqlen))
              seqid = line.strip('>\n')
              seqlen = 0
            else:
              inseq = True
              seqlen += len(line.strip(' \n'))
          
          outgff.write("##sequence-region %s 1 %d\n"%(seqid, seqlen))
          fcontig.close()
          outgff.close()
EOF
          # add the rest of the file
          tail -n +2 ${annotgff[0]}.original >> ${annotgff[0]}
        fi
        cp ${contigs}/${allcontigs} ${annotgff[0]/gff/fna}
        python ${ptgscripts}/GFFGenomeFasta2GenBankCDSProtFasta.py ${annotgff[0]} ${contigs}/${allcontigs}
      fi
      annotfna=($(ls ${annot}/${gproject}/*.fna))
      annotgbk=($(ls ${annot}/${gproject}/*.gbk))
      annotfaa=($(ls ${annot}/${gproject}/*.faa))
      annotffn=($(ls ${annot}/${gproject}/*.ffn))
      python ${ptgscripts}/add_taxid_feature2prokkaGBK.py ${annotgbk[0]} ${annotgbk[0]%*.gbk}.ptg.gbk ${straininfo}
      echo "done."
      echo "# $(dateprompt)"
    else
      echo "Error: missing mandatory information about custom genomes"
      echo "Please fill information in '${straininfo}' file before running this step again."
      exit 1
    fi
  done

  # create assembly directory and link/create relevant files
  export gblikeass=${customassemb}/genbank-format_assemblies
  mkdir -p $gblikeass/
  for gproject in `ls ${annot}` ; do
  gff=$(ls ${annot}/${gproject}/ | grep 'ptg.gff')
  assemb=${gproject}.1_${gff[0]%*.ptg.gff}
  assembpathdir=${gblikeass}/${assemb}
  assembpathrad=${assembpathdir}/${assemb}
  mkdir -p ${assembpathdir}/
  ls ${assembpathdir}/ -d
  gffgz=${assembpathrad}_genomic.gff.gz
  gzip -c ${annot}/${gproject}/$gff > $gffgz
  ls $gffgz
  gbk=$(ls ${annot}/${gproject}/ | grep 'ptg.gbk') ; gbkgz=${assembpathrad}_genomic.gbff.gz
  gzip -c ${annot}/${gproject}/$gbk > $gbkgz
  ls $gbkgz
  faa=$(ls ${annot}/${gproject}/ | grep '.faa') ; faagz=${assembpathrad}_protein.faa.gz
  gzip -c ${annot}/${gproject}/$faa > $faagz
  ffn=$(ls ${annot}/${gproject}/ | grep '.ffn') ; fnagz=${assembpathrad}_cds_from_genomic.fna.gz
  python ${ptgscripts}/format_prokkaCDS.py ${annot}/${gproject}/$ffn $fnagz
  ls $fnagz
done

ln -s ${gblikeass}/* ${indata}/assemblies/
fi

fi
#### end of the block treating custom genome set

### Groom data
## generate assembly statistics to verify genome finishing status
${ptgscripts}/groom_refseq_data.sh ${indata}/assemblies ${indata}/assembly_stats

manuin=${genomeinfo}/manual_input_metadata
mkdir -p ${manuin}/
touch ${manuin}/manual_metadata_dictionary.tab ${manuin}/manual_curated_metadata_dictionary.tab ${manuin}/manual_dbxrefs.tab
rm -rf ${genomeinfo}/assemblies*
ls ${assemblies}/* -d > ${genomeinfo}/assemblies_list
mkdir -p ${genomeinfo}/assembly_metadata
## extract assembly/sample metadata from flat files
python ${ptgscripts}/extract_metadata_from_gbff.py --assembly_folder_list=${genomeinfo}/assemblies_list --add_raw_metadata=${manuin}/manual_metadata_dictionary.tab \
--add_curated_metadata=${manuin}/manual_curated_metadata_dictionary.tab --add_dbxref=${manuin}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism" --output=${genomeinfo}/assembly_metadata

export ngenomes=$((`wc -l ${genomeinfo}/assembly_metadata/metadata.tab | cut -d' ' -f1` - 1))

mv ${envsourcescript} ${envsourcescript}0 && \
 sed -e "s#'REPLACEngenomes'#$ngenomes#" ${envsourcescript}0 > ${envsourcescript} && \
 rm ${envsourcescript}0
echo "ngenomes=$ngenomes is recorded in init file '${envsourcescript}'"

echo "Will work with a database of $ngenomes genomes"
