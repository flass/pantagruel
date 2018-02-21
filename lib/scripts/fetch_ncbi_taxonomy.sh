set -e

email=$1
DEST=$2

user='anonymous'
pswd=$email
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
source="/pub/taxonomy"
files="taxcat.tar.gz* taxcat_readme.txt taxdump.tar.gz* taxdump_readme.txt"
mkdir -p $DEST
lftp -c "open $openparam; cd $source ; mget -O ${DEST}/ $files ; quit"

# reduce taxonomic id2name reference file complexity
cd ${DEST}
for tgz in `ls *.tar.gz` ; do
    md5sum -c $tgz.md5;
    tar -xzf $tgz;
done
grep "scientific name" names.dmp | sed -e 's/\t|\t/\t/g' | cut -f1,2 > scientific_names.dmp
