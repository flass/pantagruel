#!/bin/bash

if [ -z "$3" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder 'remote_host:remote_ptg_root_folder' [files to sync ...]" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export remoteptgroot="$3"

shift 3

export ptgdb=${ptgroot}/${ptgdbname}
export envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh

remotehost=$(echo "$remoteptgroot" | cut -d':' -f1)
remotefolder=$(echo "$remoteptgroot" | cut -d':' -f2)

# create remote root folder
ssh ${remotehost} "mkdir -p ${remotefolder}"
remoteenvsourcescript=${remotefolder}/${ptgdbname}/$(basename ${envsourcescript})
# sync and translate environment source script
rsync -avz ${envsourcescript} ${remotehost}:${remoteenvsourcescript}.remote_origin
ssh ${remotehost} "sed -e 's#${ptgroot}#${remotefolder}#g' ${remoteenvsourcescript}.remote_origin ${remoteenvsourcescript}"

for sourcef in ${@} ; do
  targetf=$(echo ${sourcef%*/} | sed -e 's#${ptgroot}#${remotefolder}#')
  targetpar=$(dirname $targetf)
  ssh ${remotehost} "mkdir -p ${targetpar}"
  rsync -avz ${sourcef%*/} ${remotehost}:${targetf}

done
