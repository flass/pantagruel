#! /bin/bash
assemblies=$1
stats=$2

cd ${assemblies}/
mkdir -p ${stats}/
if [ $? != 0 ] ; then 
  echo "unable to create output directory '${stats}' ; exit now"
  exit 1
fi

filetag="_assembly_stats.txt"
patasslev="Assembly level"
patseqteq="Sequencing technology"
patcontig="all\tall\tall\tall\tcontig-"	
# assembly status
sedfilter0="s@\(.\+\)/.\+${filetag}:# ${patasslev}: \(.\+\)@\1\t\2@g"
grep "$patasslev" ./*/*${filetag} 2> /dev/null | sed -e "$sedfilter0" | sed -e 's/\r//g' > ${stats}/assembly_level.tab
# sequencing technology
sedfilter1="s@\(.\+\)/.\+${filetag}:# ${patseqteq}: \(.\+\)@\1\t\2@g"
sedfilter2="s@\(.\+\)/.\+${filetag}@\1\tNA@g"
grep "$patseqteq" ./*/*${filetag} 2> /dev/null | sed -e "$sedfilter1" | sed -e 's/\r//g' > ${stats}/sequencing_technology.tab
grep -L "$patseqteq" ./*/*${filetag} 2> /dev/null | sed -e "$sedfilter2" | sed -e 's/\r//g' >> ${stats}/sequencing_technology.tab
# contigs
for info in count N50 ; do
  sedfilter3="s@\(.\+\)/\1${filetag}:${patcontig}${info}\t\(.\+\)@\1\t\2@g"
  grep -P $patcontig ./*/*${filetag} 2> /dev/null | sed -e $sedfilter3 | grep -v $filetag > ${stats}/contig-${info}
done
