#!/bin/bash

# logging variable and function
promptdate () {
  echo $(date +'[%Y-%m-%d %H:%M:%S]') $1
}
export -f promptdate
export datepad="                      "

# function for checking that the current software version is the same as used to create the config file
checkptgversion (){
	vmsg="the current version of pantagruel (commit ${ptgversion:0:7}) is different from the one used to generate the config file '${ptgdb}/environ_pantagruel_${ptgdbname}.sh' (commit ${ptgversinit:0:7})."
  if [ "$ptgversion" != "$ptgversinit" ] ; then
    if [ "${runmode}" == 'force' ] ; then
      echo "WARNING: $vmsg"
    else
      echo "ERROR: $vmsg"
      echo "Please regenerate the config file with \`pantagruel init\` to ensure compatibility; for the same parameters to be set, just run the same command with same options as previously."
      exit 1
    fi
  fi
}

# function for checking that no pre-existing data will be over-writen, unless specified
checkfoldersafe (){
  if [ -d ${1} ] ; then
    if [ "${resumetask}" == 'true' ] ; then
      echo "Task folder '${1}' already exists; -R|--resume option was used so Pantagruel will atempt to resume from an interupted previous run"
    elif [ "${runmode}" == 'force' ] ; then
      echo "Task folder '${1}' already exists; FORCE mode is on: ERASE and recreate the folder to write new result in its place"
      rm -rf ${1}
      mkdir ${1}
    else
      echo "Task folder '${1}' already exists; will stop here rather then overwritting data."
      echo "If you want previous data to be erased, use FORCE mode with -F|--FORCE option."
      exit 2
    fi
  else
    echo "Create new task folder '${1}'"
    mkdir ${1}
  fi
}
export -f checkfoldersafe

# function for checking the success of every step
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

# function to check the presence of value in parsing option arguments
testmandatoryarg (){
  if [ -z "${2}" ]; then
   echo "ERROR: missing argument for option '${1}'" 1>&2
   echo "see pantagruel --help for more details" 1>&2
   exit 1
  fi
}


alias panup="cd ${ptgrepo} && git pull --recurse-submodules && cd -"