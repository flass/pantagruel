#!/bin/bash

assembly_folder=$1
DEST=$2

cd ${assembly_folder}
mkdir -p ${DEST}
filetag="_assembly_stats.txt"
patasslev="Assembly level"
patseqteq="Sequencing technology"
patcontig="all\tall\tall\tall\tcontig-"	
# assembly status
sedfilter0="s@\(.\+\)/.\+${filetag}:# ${patasslev}: \(.\+\)@\1\t\2@g"
grep "$patasslev" */*${filetag} | sed -e "$sedfilter0" | sed -e 's/\r//g' > ${DEST}/assembly_level.tab
# assembly status
sedfilter0="s@\(.\+\)/.\+${filetag}:# ${patasslev}: \(.\+\)@\1\t\2@g"
grep "$patasslev" */*${filetag} | sed -e "$sedfilter0" | sed -e 's/\r//g' > ${DEST}/assembly_level.tab
# sequencing technology
sedfilter1="s@\(.\+\)/.\+${filetag}:# ${patseqteq}: \(.\+\)@\1\t\2@g"
sedfilter2="s@\(.\+\)/.\+${filetag}@\1\tNA@g"
grep "$patseqteq" */*${filetag} | sed -e "$sedfilter1" | sed -e 's/\r//g' > ${DEST}/sequencing_technology.tab
grep -L "$patseqteq" */*${filetag} | sed -e "$sedfilter2" | sed -e 's/\r//g' >> ${DEST}/sequencing_technology.tab
# contigs
for info in count N50 ; do
  sedfilter3="s@\(.\+\)/\1${filetag}:${patcontig}${info}\t\(.\+\)@\1\t\2@g"
  grep -P $patcontig */*${filetag} | sed -e $sedfilter3 | grep -v $filetag > ${DEST}/contig-${info}
done
exit 0
