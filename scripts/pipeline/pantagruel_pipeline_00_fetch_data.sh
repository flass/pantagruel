#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

parseGBass (){
python2.7 << EOF
import re
gbassdir = '${1}'
assembpatgenbank = re.compile('(GC[AF]_[^\._]+\.[0-9])_(.+)')
assembsearch = assembpatgenbank.search(gbassdir)
if assembsearch:
  assacc, assname = assembsearch.groups()
  print "%s\t%s"%(assacc, assname)

EOF
}
export -f parseGBass

extractass (){
  srcass=${1}
  lndestass=${2}
  gp2ass=${3}
  ## extract content of a folder of assemblies as those obtained using NCBI web interface:
  # Search Assembly database using a query defined as convenient:  
  # e.g.: 'txid543[Organism:exp] AND ("latest refseq"[filter] AND all[filter] NOT anomalous[filter]) AND ("complete genome"[filter])'
  # save assemblies (Download Assemblies > Source Database: RefSeq; File type: All file type (including assembly-structure directory))
  # as ${srcass}/genome_assemblies.tar or as ${srcass}/genome_assemblies_DATASETTAG.tar if several datasets have to be merged
  echo "extract assembly data from folder '${srcass}'"
  cd ${srcass}/
  # look for assembly folders
  ls -A ${srcass}/ | grep -v 'assemblies' | grep -v 'md5checksum' > ${srcass}/genome_assemblies_list
  # look for archives containing many assembly folders (datasets downloaded from NCBI Assembly website)
  assarch=$(ls genome_assemblies*.tar 2> /dev/null)
  if [ "$assarch" ] ; then
    echo "detected archive(s) of genome assemblies (presumably downloaded from NCBI Assembly website):"
    echo "$assarch"
    # extract assembly folders from the archives
    for tarf in `ls genome_assemblies*.tar` ; do 
      tartag=${tarf#genome_assemblies*}
      datasettag=${tartag%.*}
      tar -tf genome_assemblies${datasettag}.tar | cut -d'/' -f2 | sort -u | grep -v "README\|report" > ${srcass}/genome_assemblies${datasettag}_list
      tar -xf genome_assemblies${datasettag}.tar 
      [ ! -z "${datasettag}" ] && mv -f report.txt report${datasettag}.txt
      assd=(`ls ncbi-genomes-* -d`)
      for dass in `ls ${assd[0]}/ | grep -v README` ; do
        if [ ! -z "${dass}" ] ; then
          rm -rf ${srcass}/$dass
          mv ${assd[0]}/$dass ${srcass}/
		  if [ ! -z "${gp2ass}" ] ; then
		     assaccname="$(parseGBass ${dass})"
			 if [ -z "${assaccname}" ] ; then 
			   echo "could not parse '${dass}' assembly name as expected from the NCBI Assembly format pattern (Accession.v_Name => python regex: '(GC[AF]_[^\._]+\.[1-9])_(.+)')"
			 fi
	         echo -e "${dass}\t${assaccname}\tNCBI_Assembly_extracted" >> ${gp2ass}
		  fi
        fi
      done
      rm -r ${assd[0]}/ 
    done
  fi
  grep "# Organism name" ${srcass}/*/*_assembly_stats.txt > ${srcass}/all_assemblies_organism_names
  
  # record links in centralised folder of assemblies
  mkdir -p ${lndestass}/
  #~ relpathass2srcass=$(realpath --relative-to=${lndestass} ${srcass})
  #~ for ass in `cat ${srcass}/genome_assemblies*_list` ; do
    #~ ln -s ${relpathass2srcass}/${ass} ${lndestass}/
  #~ done
  # better preserve actual path as input data are supposed to be located outside the Pantagruel database
  for ass in `cat ${srcass}/genome_assemblies*_list` ; do
    ln -s ${srcass}/${ass} ${lndestass}/
  done
  cd -
}
export -f extractass

downloadass (){
  srclist=${1}
  lndestass=${2}
  dldestass=${3}
  gp2ass=${4}
  echo "$(promptdate) fetch assembly data from NCBI FTP accordng to list '${srclist}'"
  # fetch genome assemblies from NCBI FTP based on list prvided with -L option
  if [ ! -z "${dldestass}" ] ; then
    srcassftpdest=${dldestass}/$(basename ${srclist})_assemblies_from_ftp
  else
    srcassftpdest=${srclist}_assemblies_from_ftp
  fi
  mkdir -p ${srcassftpdest}/
  ${ptgscripts}/fetch_refseq_assemblies_from_list.sh ${srclist} ${srcassftpdest}
  checkexec "could not fetch all the genomes ; exit now"
  
  # record links in centralised folder of assemblies
  mkdir -p ${lndestass}/
  # better link with relative path as input data are suposed to be dowloaded inside the Pantagruel database
  relpathass2srcass=$(realpath --relative-to=${lndestass} ${srcassftpdest})
  for ass in `ls ${srcassftpdest}/` ; do
    ln -s ${relpathass2srcass}/${ass} ${lndestass}/
    if [ ! -z "${gp2ass}" ] ; then
		 assaccname=$(parseGBass ${ass})
		 if [ -z "${assaccname}" ] ; then 
		   echo "could not parse '${ass}' assembly name as expected from the NCBI Assembly format pattern (Accession.v_Name => python regex: '(GC[AF]_[^\._]+\.[1-9])_(.+)')"
		 fi
		 echo -e "${ass}\t${assaccname}\tNCBI_Assembly_downloaded" >> ${gp2ass}
    fi
  done
}
export -f downloadass

parsefastaext (){
python2.7 << EOF
allcontigs = '${1}'
for ext in ['.fa', '.fna', '.fasta', '.fsa']:
  if allcontigs.endswith(ext):
    print allcontigs.rsplit(ext, 1)[0]
    break
else:
  print allcontigs

EOF
}
export -f parsefastaext

addregionannotfeat (){
python2.7 << EOF
fcontig = open('${1}', 'r')
outgff = open('${2}', 'a')
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
}
export -f addregionannotfeat

checkptgversion
checkfoldersafe ${indata}

#################################
## 00. Data download and grooming
#################################

if [ ! -e ${ncbitax}/names.dmp ] ; then
echo "$(promptdate) did not find the relevant NCBI Taxonomy flat files in '${ncbitax}/'; download the from NCBI Taxonomy FTP"
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
  #~ lftp ${openparam} -e "cd ${source} ; mget -O ${ncbitax}/ ${files} ; quit"
  lftp -c "set dns:fatal-timeout 10 ; open ${openparam} ; cd ${source} ; mget -O ${ncbitax}/ ${files} ; quit"
  checkexec "could not download the NCBI Taxonomy files; please check your network connection or download manually the files '${files}' from the FTP site ftp://${host}${source} into folder '${ncbitax}/' ; exit now" "succesfully downloaded  NCBI Taxonomy flat files in '${ncbitax}/'"
  cd ${ncbitax}/
  for tgz in `ls *.tar.gz` ; do md5sum -c ${tgz}.md5 && tar -xzf ${tgz} ; done
  cd -
fi

mkdir -p ${assemblies}/
echo -e "genome_source\tassembly_id\tassembly_name\torigin" > ${gp2ass}

if [ ! -z "${ncbiass}" ] ; then
  extractass ${ncbiass} ${assemblies} ${gp2ass}
fi

if [ ! -z "${listncbiass}" ] ; then
  downloadass ${listncbiass} ${assemblies} ${indata} ${gp2ass}
fi

if [ ! -z "${refass}" ] ; then
  extractass ${refass} ${prokkaref}
fi

if [ ! -z "${listrefass}" ] ; then
  downloadass ${listrefass} ${prokkaref} ${indata}
fi


if [ ! -z "${customassemb}" ] ; then
  echo "$(promptdate) extract assembly data from folder '${customassemb}'"
  
  contigsets="$(ls -A "${contigs}/" 2>/dev/null)"
  echo "found $(echo $contigsets | wc -w) contig files (raw genome assemblies) in ${contigs}/"
  ## if the folder of custom/user-provided set of genomes is not empty
  if [[ -z "${contigsets}" ]] ; then
    echo "folder ${contigs} is empty! did you intend to provide custom genomeassemblies with -a option? ; exit now"
    exit 1
  else
    
    # gather representative proteins from the custom reference genome set to make a prot database for Prokka to search for similarities
    prevrefdb="$(${ptgscripts}/make_prokka_ref_genus_db.sh 'check' ${refgenus})"
    if [[ "${resumetask}" == 'true' && ! -z "${prevrefdb}" ]] ; then
      echo "Resume mode: the reference database for ${refgenus} already exsists:"
      ls ${prevrefdb}
      echo "skip the reference db building"
    else
      if [ ! -z "$(ls -A "${prokkaref}/" 2>/dev/null)" ] ; then
        echo "$(promptdate) generate Prokka reference database for annotation of genomes from the genus '${ptgscripts}/make_prokka_ref_genus_db.sh ${prokkaref} ${refgenus} ${ptgtmp} ${ptglogs}' based on assemblies found in '${prokkaref}' (specified through options --refseq_ass4annot or --refseq_list4annot)"
        ${ptgscripts}/make_prokka_ref_genus_db.sh ${prokkaref} ${refgenus} ${ptgtmp} ${ptglogs}
      elif [ ! -z "$(ls -A "${assemblies}/" 2>/dev/null)" ] ; then
        echo "$(promptdate) generate Prokka reference database for annotation of genomes from the genus '${refgenus}' based on assemblies found in '${assemblies}' (specified through options -A or -L)"
        ${ptgscripts}/make_prokka_ref_genus_db.sh ${assemblies} ${refgenus} ${ptgtmp} ${ptglogs}
      fi
    fi
  
    mkdir -p ${annot}
	mkdir -p ${gblikeass}/
	
    for allcontigs in ${contigsets} ; do
#      gproject=${allcontigs%%.fa*}
	  gproject=$(parsefastaext ${allcontigs})
      echo "$(promptdate) ${gproject}"
      if [ -d ${custannot}/${gproject} ] ; then
        echo "found annotation folder '${custannot}/${gproject}' ; skip annotation of contigs in '${contigs}/${allcontigs}'"
        # better preserve actual path as input data are external to the Pantagruel database 
        ln -s ${custannot}/${gproject} ${annot}/
      else
	    if [[ "${resumetask}" == 'true' && -d ${annot}/${gproject} && ! -z $(ls ${annot}/${gproject}/*.gff) && ! -z $(ls ${annot}/${gproject}/*.gbk) ]] ; then
		  echo "found already computed annotation in ${annot}/${gproject}/:"
		  rm -f ${annot}/${gproject}/*.ptg.gbk ${annot}/${gproject}/*.ptg.gff
		  ls -ltr ${annot}/${gproject}/
		  echo "skip running Prokka"
	    elif [[ "${resumetask}" == 'true' && -s ${annot}/${gproject}.tar.gz && "$(tar -tf ${annot}/${gproject}.tar.gz | grep -v '\.ptg\.' | grep -c '\.gff\|\.gbk')" -eq 2 ]] ; then
		  echo "found already computed annotation in archive ${annot}/${gproject}.tar.gz that already contains .gff and .gbk files"
		  echo "skip running Prokka"
		else  
          echo "will annotate contigs in '${contigs}/${allcontigs}'"
          promptdate
          # run Prokka 
          echo "### assembly: $gproject; contig files from: ${contigs}/${allcontigs}"
          echo "running Prokka..."
          ${ptgscripts}/run_prokka.sh ${gproject} ${contigs}/${allcontigs} ${straininfo} ${annot}/${gproject} ${seqcentre} &> ${ptglogs}/${ptgdbname}_customassembly_annot_prokka.${gproject}.log
          checkexec "something went wrong when annotating genome '${gproject}'; check log at '${ptglogs}/${ptgdbname}_customassembly_annot_prokka.${gproject}.log'" "succesfully annotated genome '${gproject}'"
          promptdate
        fi
	  fi
	  doprocess1=true
      if [[ "${resumetask}" == 'true' && ! -d ${annot}/${gproject} && -s ${annot}/${gproject}.tar.gz ]] ; then
	    if [ "$(tar -tf ${annot}/${gproject}.tar.gz | grep -c '.ptg.gff\|.ptg.gbk')" -eq 2 ] ; then
	      echo "annotation in archive ${annot}/${gproject}.tar.gz already contains .ptg.gff and .ptg.gbk files"
	     doprocess1=false
		else
		  echo "annotation in archive ${annot}/${gproject}.tar.gz does not contain .ptg.gff and .ptg.gbk files; extract fill to process the input .gff and .gbk files"
		  tar -xzf ${annot}/${gproject}.tar.gz
		fi
	  fi
	  if [[ "${doprocess1}" == 'true' ]] ; then
        if [[ "$(grep ${gproject} ${straininfo} | cut -f2- | grep -P '[^\t]' | wc -l)" -eq 0 ]] ; then
          echo "Error: missing mandatory information about custom genomes"
          echo "Please fill information in '${straininfo}' file before running this step again."
          exit 1
        fi
        echo "fix annotation to integrate region information into GFF files"
        annotgff=($(ls ${annot}/${gproject}/*.gff))
        if [ -z ${annotgff[0]} ] ; then
          echo "Error: missing mandatory GFF file in ${annot}/${gproject}/ folder"
          exit 1
        fi
        if [ -z "$(grep '##sequence-region' ${annotgff[0]})" ] ; then
        mv -f ${annotgff[0]} ${annotgff[0]}.original
          # make GFF source file look more like the output of Prokka
          head -n1 ${annotgff[0]}.original > ${annotgff[0]}
          # add region annotation features
          addregionannotfeat ${contigs}/${allcontigs} ${annotgff[0]}
          # add the rest of the file, in order of contig names!
          tail -n +2 ${annotgff[0]}.original | sort >> ${annotgff[0]}
        fi
        mkdir -p ${ptglogs}/add_region_feature2prokkaGFF/
		addregfeatlog=${ptglogs}/add_region_feature2prokkaGFF/add_region_feature2prokkaGFF.${gproject}.log
		annotptggff=${annotgff[0]/.gff/.ptg.gff}
        python2.7 ${ptgscripts}/add_region_feature2prokkaGFF.py ${annotgff[0]} ${annotptggff} ${straininfo} ${contigs}/${allcontigs} ${assembler} &> ${addregfeatlog}
        checkexec "something went wrong when adding region features to GFF file (In: ${annotgff[0]}; Out:${annotptggff}; Stdin/Stderr or the last command: ${addregfeatlog})"
        echo "fix annotation to integrate taxid information into GBK files"
        annotfna=($(ls ${annot}/${gproject}/*.fna))
        annotgbk=($(ls ${annot}/${gproject}/*.gbf 2> /dev/null || ls ${annot}/${gproject}/*.gbk))
        annotfaa=($(ls ${annot}/${gproject}/*.faa))
        annotffn=($(ls ${annot}/${gproject}/*.ffn))
        if [[ -z "${annotfna[0]}" || -z "${annotgbk[0]}" || -z "${annotffn[0]}" || -z "${annotfaa[0]}" ]] ; then
          echo "At least one of these files is missing in ${annot}/${gproject}/ folder: contig fasta file (.fna), GenBank flat file (gbk/gbf), CDS Fasta (ffn) or protein Fasta (faa)."
          echo "Will (re)generate them from the GFF anotation and genomic Fasta sequence; files already present are kept with an added prefix '.original'"
          for annotf in ${annotfna[@]} ${annotgbk[@]} ${annotffn[@]} ${annotfaa[@]}; do
            if [[ ! -z "${annotf}" ]] ; then
              mv ${annotf} ${annotf}.original
            fi
          done
          cp ${contigs}/${allcontigs} ${annotgff[0]/gff/fna}
          python2.7 ${ptgscripts}/GFFGenomeFasta2GenBankCDSProtFasta.py ${annotptggff} ${annotgff[0]/gff/fna}
        checkexec "something went wrong when generating the GenBank flat file from GFF file ${annotptggff}" "succesfuly generated the GenBank flat file from GFF file ${annotptggff}"
        fi
        annotfna=($(ls ${annot}/${gproject}/*.fna))
        annotgbk=($(ls ${annot}/${gproject}/*.gbk 2> /dev/null))
        annotrad=${annotgbk[0]%*.gbk}
        if [ -z "${annotgbk[0]}" ] ; then
          annotgbk=($(ls ${annot}/${gproject}/*.gbf))
          annotrad=${annotgbk[0]%*.gbf}
        fi
        annotfaa=($(ls ${annot}/${gproject}/*.faa))
        annotffn=($(ls ${annot}/${gproject}/*.ffn))
        python2.7 ${ptgscripts}/add_taxid_feature2prokkaGBK.py ${annotgbk[0]} ${annotrad}.ptg.gbk ${straininfo}
        checkexec "something went wrong when modifying the GenBank flat file ${annotgbk[0]}" "succesfully modified the GenBank flat file ${annotgbk[0]}"
	  fi
#    done
  
      # create assembly directory and link/create relevant files
      echo "will create GenBank-like assembly folder for user-provided genomes"
#    mkdir -p ${gblikeass}/
#    for gproject in `ls ${annot}` ; do
	  doprocess2=true
      gff=$(ls ${annot}/${gproject}/ 2> /dev/null | grep 'ptg.gff' | grep -v '\.original')
	  if [ -z "${gff}" ] ; then
	    if [ -s ${annot}/${gproject}.tar.gz  ] ; then
		  # annotation was previously processed then compressed; try and see if the final files are present
		  gblikefilemissing=-1
		  assembpathdir=$(ls -d ${gblikeass}/${gproject}.1_*)
		  if [ ! -z "${assembpathdir}" ] ; then
		    gblikefilemissing=$(( ${gblikefilemissing} + 1 ))
	        for gbext in '_cds_from_genomic.fna.gz' '_genomic.fna.gz' '_genomic.gbff.gz' '_genomic.gff.gz' '_protein.faa.gz' ; do
			  if [ -z "$(ls -d ${gblikeass}/${gproject}.1_*/*${gbext})" ] ; then
			    gblikefilemissing=$(( ${gblikefilemissing} + 1 ))
		      fi
		    done
		  fi
		  if [ ${gblikefilemissing} -eq 0 ] ; then
		    echo "all final GenBank-like files found in folder ${assembpathdir}/ ; skip processing."
			doprocess2=false
			assemb=$(basename ${assembpathdir})
		  else
		    tar -xzf ${annot}/${gproject}.tar.gz && gff=$(ls ${annot}/${gproject}/ 2> /dev/null | grep 'ptg.gff' | grep -v '\.original')
	        if [ -z "${gff}" ] ; then
	          echo "Error: could not find file '${annotptggff}' in the extracted content of archive '${annot}/${gproject}.tar.gz'; annotation is missing for genome '${gproject}': exit now"
		      exit 1
			else
			  assemb=${gproject}.1_${gff[0]%*.ptg.gff}
	        fi
	      fi
		else
		  echo "Error: could not find file '${annotptggff}' or archive '${annot}/${gproject}.tar.gz' ; annotation is missing for genome '${gproject}': exit now"
		  exit 1
	    fi
	  else
	    assemb=${gproject}.1_${gff[0]%*.ptg.gff}
	  fi
      
	  assaccname=$(parseGBass ${gproject})
	  if [ ! -z "${assaccname}" ] ; then
	    echo -e "${gproject}\t${assaccname}.1_${gff[0]%*.ptg.gff}\tcustom_assembly" >> ${gp2ass}
      else
	    echo -e "${gproject}\t${gproject}.1\t${gff[0]%*.ptg.gff}\tcustom_assembly" >> ${gp2ass}
	  fi
	  if [[ "${doprocess2}" == 'true' ]] ; then
        assembpathdir=${gblikeass}/${assemb}
        assembpathrad=${assembpathdir}/${assemb}
        mkdir -p ${assembpathdir}/
        ls -l ${assembpathdir}/ -d
        gffgz=${assembpathrad}_genomic.gff.gz
        gzip -c ${annot}/${gproject}/${gff} > ${gffgz}
        ls -l ${gffgz}
        gbk=($(ls ${annot}/${gproject}/ | grep 'ptg.gbk' | grep -v '.original'))
        gbkgz=${assembpathrad}_genomic.gbff.gz
        gzip -c ${annot}/${gproject}/${gbk[0]} > ${gbkgz}
        ls -l ${gbkgz}
        faa=($(ls ${annot}/${gproject}/ | grep '.faa' | grep -v '.original'))
        faagz=${assembpathrad}_protein.faa.gz
        gzip -c ${annot}/${gproject}/${faa[0]} > ${faagz}
        ls -l ${faagz}
        gfna=($(ls ${annot}/${gproject}/ | grep '.fna' | grep -v '.original'))
        gfnagz=${assembpathrad}_genomic.fna.gz
        gzip -c ${annot}/${gproject}/${gfna[0]} > ${gfnagz}
        ls -l ${gfnagz}
        ffn=($(ls ${annot}/${gproject}/ | grep '.ffn' | grep -v '.original'))
        cdsfnagz=${assembpathrad}_cds_from_genomic.fna.gz
        python2.7 ${ptgscripts}/format_prokkaCDS.py ${annot}/${gproject}/${ffn[0]} ${cdsfnagz}
        ls -l ${cdsfnagz}
	    if [[ "${compress}" == 'on' ]] ; then
            # compress and only upon success delete the uncompress folder
		echo "will compress folder '${annot}/${gproject}/':"
          tar -czf ${annot}/${gproject}.tar.gz ${annot}/${gproject}/ && rm -rf ${annot}/${gproject}/ && echo -e "successfully compressed into\n:'$(ls -l ${annot}/${gproject}.tar.gz)' and deleted source folder." || "Compression failed; leave source folder '${annot}/${gproject}/' as is"
	    fi
	  fi
#    done > ${ptglogs}/genbank-format_assemblies.log
    done
    
    relpathass2gblass=$(realpath --relative-to=${assemblies} ${gblikeass})
    for gblass in $(ls -A ${gblikeass}/) ; do
      ln -s ${relpathass2gblass}/${gblass} ${assemblies}/
    done
  fi
fi

#### end of the block treating custom genome set

#### Summarize data
mkdir -p ${genomeinfo}/
rm -rf ${genomeinfo}/assemblies*
ls -Ad ${assemblies}/* > ${genomeinfo}/assemblies_list

### Groom data
## test that there are at least 4 genomes to process
export ngenomes=$(ls -A "${assemblies}/" 2>/dev/null | wc -l)
if [ ${ngenomes} -lt 4 ] ; then
  echo "Error: there are not enough genomes (< 4) in input for pantegruel to produce anything meaningful."
  echo "please add more genome assemblies in either ${ncbiass}/ or ${contigs}/ (see manual with \`pantagruel -- help\`)"
  echo "exit now."
  exit 1
fi

## verify the presence of '*_cds_from_genomic.fna[.gz]' file, as it is sometimes missing even from a RefSeq-downloaded assembly folder (happens for recently published assemblies)
mkdir -p ${indata}/extracted_cds_from_genomic_fasta 
python2.7 ${ptgscripts}/check_create_cds_from_genomic.py ${genomeinfo}/assemblies_list ${indata}/extracted_cds_from_genomic_fasta
for dass in $(ls -A ${indata}/extracted_cds_from_genomic_fasta) ; do
  relpathass2extcds=$(realpath --relative-to=${assemblies}/${dass} ${indata}/extracted_cds_from_genomic_fasta/${dass})
  for nfcds in $(ls -A ${indata}/extracted_cds_from_genomic_fasta/${dass}/${dass}_cds_from_genomic.fna*) ; do
    bnnfcds=$(basename ${nfcds})
    ln -s ${relpathass2extcds}/${bnnfcds} ${assemblies}/${dass}/
  done
done

## generate assembly statistics to verify genome finishing status
${ptgscripts}/groom_refseq_data.sh ${assemblies} ${indata}/assembly_stats

manuin=${genomeinfo}/manual_input_metadata
mkdir -p ${manuin}/
touch ${manuin}/manual_metadata_dictionary.tab ${manuin}/manual_curated_metadata_dictionary.tab ${manuin}/manual_dbxrefs.tab
mkdir -p ${genomeinfo}/assembly_metadata
## extract assembly/sample metadata from flat files
python2.7 ${ptgscripts}/extract_metadata_from_gbff.py --assembly_folder_list=${genomeinfo}/assemblies_list --add_raw_metadata=${manuin}/manual_metadata_dictionary.tab \
--add_curated_metadata=${manuin}/manual_curated_metadata_dictionary.tab --add_dbxref=${manuin}/manual_dbxrefs.tab --add_assembly_info_dir=${indata}/assembly_stats \
--default_species_name="unclassified organism" --output=${genomeinfo}/assembly_metadata


## compute genome-to-genome MASH distances and plot them, notably as a heatmap along a distance tree (and possibly along the core-genome reference tree)
if [ ! -z "$(command -v mash)" ] ; then 
  if [ ! -z "${ptgthreads}" ] ; then 
    paramash="-p ${ptgthreads}"
  else
    paramash="-p $(nproc)"
  fi
  fnalist=${indata}/all_assemblies_genomic_fasta_list
  rm -f ${fnalist}
  for ass in $(ls ${assemblies}) ; do
    ls ${assemblies}/${ass}/${ass}_genomic.fna.gz 2> /dev/null >> ${fnalist} || ls ${assemblies}/${ass}/${ass}_genomic.fna 2> /dev/null >> ${fnalist} 
  done
  mash triangle ${paramash} $(cat ${fnalist}) > ${indata}/all_assemblies_mash.dist
  
  if [[ -s "${speciestree}" && -s "${database}/genome_codes.tab" ]] ; then
    # only likely to happen if task 00 is re-run with -R after tasks 03 and 05 are complete
    Rscript ${ptgscripts}/plotmashdistcluster.r ${indata}/all_assemblies_mash.dist ${genomeinfo}/assembly_metadata/metadata.tab ${speciestree} ${database}/genome_codes.tab
  else
    Rscript ${ptgscripts}/plotmashdistcluster.r ${indata}/all_assemblies_mash.dist ${genomeinfo}/assembly_metadata/metadata.tab
  fi
  
fi
promptdate
