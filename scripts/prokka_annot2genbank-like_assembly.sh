#!/usr/bin/env bash

# define functions
promptdate () {
  echo $(date +'[%Y-%m-%d %H:%M:%S]') $1
}
export -f promptdate
export datepad="                      "

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

checkexec (){
  if [ ${?} != 0 ]; then
    echo "ERROR: $1" 1>&2
    exit 1
  else
    if [ ! -z "$2" ] ; then
      echo -e "$2"
    fi
  fi
}
export -f checkexec
 
#parseGBass (){
#python2.7 << EOF
#import re
#gbassdir = '${1}'
#assembpatgenbank = re.compile('(GC[AF]_[^\._]+\.[0-9])_(.+)')
#assembsearch = assembpatgenbank.search(gbassdir)
#if assembsearch:
#  assacc, assname = assembsearch.groups()
#  print "%s\t%s"%(assacc, assname)
#
#EOF
#}
#export -f parseGBass

# parse script arguments
[ -z "${3}" ] && echo "Error: missing arguments" && echo -e "usage is: ${0} contig_folder_(original_prokka_input) input_prokka_annotation_folder output_GenBank-like_annotation_folder strain_info_file [pantagruel_scripts_folder] [log_folder]"
contigs=${1}
annot=${2}
gblikeass=${3}
straininfo=${4}
[ -z "${ptgscripts}" ] && ptgscripts=${5}
[ -z "${ptgscripts}" ] && echo "Error: please define \${ptgscripts} environment variable or pass its value as 5th positional argument ; exit now" && exit 1
[ -z "${ptglogs}" ] && ptglogs=${PWD}/${0/.sh/_logs}

mkdir -p ${gblikeass}/
mkdir -p ${ptglogs}/
contigsets="$(ls -A "${contigs}/" 2>/dev/null)"

nannnot=0

for allcontigs in ${contigsets} ; do
  gproject=$(parsefastaext ${allcontigs})
  echo "$(promptdate) ${gproject}"
  doprocess1=true
  if [[ ! -d ${annot}/${gproject} && -s ${annot}/${gproject}.tar.gz ]] ; then
	if [ "$(tar -tf ${annot}/${gproject}.tar.gz | grep -c '.ptg.gff\|.ptg.gbk')" -eq 2 ] ; then
	  echo "annotation in archive ${annot}/${gproject}.tar.gz already contains .ptg.gff and .ptg.gbk files"
	 doprocess1=false
	else
	  echo "annotation in archive ${annot}/${gproject}.tar.gz does not contain .ptg.gff and .ptg.gbk files; extract fill to process the input .gff and .gbk files"
	  tar -xzf ${annot}/${gproject}.tar.gz
	fi
  elif [[ ! -d ${annot}/${gproject} ]] ; then
    echo "Error: annotation folder or archive is missing for contig set '${allcontigs}' ; exit now" && exit 1
  fi
  if [[ "${doprocess1}" == 'true' ]] ; then
	if [[ "$(grep ${gproject} ${straininfo} | cut -f2- | grep -P '[^\t]' | wc -l)" -eq 0 ]] ; then
	  echo "Error: missing mandatory information about custom genome assembly with id '${gproject}'"
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
	python2.7 ${ptgscripts}/add_region_feature2prokkaGFF.py ${gproject} ${annotgff[0]} ${annotptggff} ${straininfo} ${contigs}/${allcontigs} ${assembler} &> ${addregfeatlog}
	checkexec "something went wrong when adding region features to GFF file (In: ${annotgff[0]}; Out:${annotptggff}; Stdin/Stderr or the last command: ${addregfeatlog})"
	annotfna=($(ls ${annot}/${gproject}/*.fna))
	annotgbk=($(ls ${annot}/${gproject}/*.gbf 2> /dev/null || ls ${annot}/${gproject}/*.gbk))
	annotfaa=($(ls ${annot}/${gproject}/*.faa))
	annotffn=($(ls ${annot}/${gproject}/*.ffn))
	if [[ -z "${annotfna[0]}" || -z "${annotgbk[0]}" || -z "${annotffn[0]}" || -z "${annotfaa[0]}" ]] ; then
	  echo "Error: At least one of these files is missing in ${annot}/${gproject}/ folder: contig fasta file (.fna), GenBank flat file (gbk/gbf), CDS Fasta (ffn) or protein Fasta (faa)."
	  exit 1
#	  echo "At least one of these files is missing in ${annot}/${gproject}/ folder: contig fasta file (.fna), GenBank flat file (gbk/gbf), CDS Fasta (ffn) or protein Fasta (faa)."
#	  echo "Will (re)generate them from the GFF anotation and genomic Fasta sequence; files already present are kept with an added prefix '.original'"
#	  for annotf in ${annotfna[@]} ${annotgbk[@]} ${annotffn[@]} ${annotfaa[@]} ; do
#		if [[ ! -z "${annotf}" ]] ; then
#		  mv ${annotf} ${annotf}.original
#		fi
#	  done
#	  cp ${contigs}/${allcontigs} ${annotgff[0]/gff/fna}
#	  python2.7 ${ptgscripts}/GFFGenomeFasta2GenBankCDSProtFasta.py ${annotptggff} ${annotgff[0]/gff/fna}
#	  checkexec "something went wrong when generating the GenBank flat file from GFF file ${annotptggff}" "succesfuly generated the GenBank flat file from GFF file ${annotptggff}"
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
	echo "fix annotation to integrate taxid information into GBK files"
	python2.7 ${ptgscripts}/add_taxid_feature2prokkaGBK.py ${gproject} ${annotgbk[0]} ${annotrad}.ptg.gbk ${straininfo}
	checkexec "something went wrong when modifying the GenBank flat file ${annotgbk[0]}" "succesfully modified the GenBank flat file ${annotgbk[0]}"
  fi

  # create assembly directory and link/create relevant files
  echo "will create GenBank-like assembly folder for user-provided genomes"
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
  nannot=$(( ${nannot} + 1 ))
done
    
echo "$(promptdate) converted ${nannot} from Prokka output format from '${annot}/' to a GenBank Assembly style format, wirtten in '${gblikeass}/'; done."