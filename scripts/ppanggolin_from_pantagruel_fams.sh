#!/usr/bin/env bash

usage (){
  echo "ppangolin_from_pantagruel_fams.sh ptgroot ptgdbname [pangotag] [ncpus]"
  echo "  ptgroot	folder location of source Pantagruel database"
  echo "  ptgdbname	name of folder containing the source Pantagruel database"
  echo "  pangotag	tag for PPaNGOLIN result folder within the Pantagruel database (default=today's date in YYYY-MM-DD format)"
  echo "  ncpus 	number of parallel cores (default=4)"
}

if [ ! -z "${1}" ] ; then
  export ptgroot=${1}
elif [ ! -z "${ptgroot}" ] ; then
  echo "using environment variable ptgroot=${ptgroot}"
else
  echo "Error: need to set environment variable \$ptgroot or pass value through \$1 argument"
  usage
  exit 1
fi
if [ ! -z "${2}" ] ; then
  export ptgdbname=${2}
elif [ ! -z "${ptgdbname}" ] ; then
  echo "using environment variable ptgdbname=${ptgdbname}"
else
  echo "Error: need to set environment variable \$ptgdbname or pass value through \$2 argument"
  usage
  exit 1
fi
if [ ! -z "${3}" ] ; then
  export pangotag=${3}
elif [ ! -z "${pangotag}" ] ; then
  echo "using environment variable pangotag=${pangotag}"
else
  export pangotag=$(date "+%y-%m-%d")
  echo "setting the PPANGOLIN database tag to default value, today's date: ${pangotag}"
fi
if [ ! -z "${4}" ] ; then
  export ncpus=${4}
elif [ ! -z "${ncpus}" ] ; then
  echo "using environment variable ncpus=${ncpus}"
else
  export ncpus=4
  echo "setting the number of parallel cores to default value: ${ncpus}"
fi

export ptgdb=${ptgroot}/${ptgdbname}

source ${ptgdb}/environ_pantagruel_${ptgdbname}.sh

export ppanggo=${ptgdb}/ppanggolin_${pangotag}
export ppangglog=${ppanggo}/logs
export ppanggoin=${ppanggo}/input
export ppanggout=${ppanggo}/pangenome
mkdir -p ${ppangglog}/ ${ppanggoin}/

# load annotated genomes
sqlite3 ${sqldb} "SELECT code, assembly_id, assembly_name from assemblies;" | sed -e 's/|/\t/g' | while read code assid assname ; do
echo -e "${code}\t$(ls ${indata}/genbank-format_assemblies/${assid}_${assname}/*gbff.gz)"
done > ${ppanggo}/input/genome_annotations.txt

ppanggolin annotate --cpu ${ncpus} --anno ${ppanggo}/input/genome_annotations.txt --output ${ppanggout} --log ${ppangglog}/annotate.log

# load gene families
sqlite3 ${sqldb} "SELECT gene_family_id, locus_tag from coding_sequences;" | sed -e 's/|/\t/g' > ${ppanggo}/input/gene_families.txt
grep -v 'VCHOLC000000' ${ppanggo}/input/gene_families.txt > ${ppanggo}/input/gene_families_noORFans.txt

ppanggolin cluster --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --clusters ${ppanggo}/input/gene_families_noORFans.txt --infer_singletons --log ${ppangglog}/cluster.log

# run the rest of the workflow
ppanggolin graph --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --remove_high_copy_number 5 --log ${ppangglog}/graph.log

ppanggolin partition --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --log ${ppangglog}/partition.log

ppanggolin rgp --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --log ${ppangglog}/rgp.log
ppanggolin write --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --regions --output ${ppanggo}/rgp

ppanggolin spot --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --log ${ppangglog}/spot.log
ppanggolin write --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --spots --output ${ppanggo}/spots

ppanggolin module --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --log ${ppangglog}/module.log
ppanggolin write --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --modules --output ${ppanggo}/modules
ppanggolin write --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --spot_modules --output ${ppanggo}/spot_modules

ppanggolin write --cpu ${ncpus} -p ${ppanggout}/pangenome.h5 --json --gexf --light_gexf --output ${ppanggo}/graph