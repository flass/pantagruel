#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

CRANmirror="https://cran.ma.imperial.ac.uk/" # should be edited to match nearest CRAN mirror

checkexec (){
  if [ $? != 0 ]; then
    echo "ERROR: $1" 1>&2
    echo "ERROR: abort installation of Pantagruel"
    exit 1
  else
    if [ ! -z "$2" ] ; then
      echo "$2"
    fi
  fi
}

usage (){
	echo "Usage: $0 software_dir [bin_dir] [--no-interpro]"
	echo "  software_dir: target installation folder for pantagruel and other software"
	echo "  bin_dir: folder where to link relevant executable files, including the main `pantagruel` excutable (defaults to ~/bin/); the path to this folder will be permatnently added to your path (via editing your ~/.bashrc)"
	echo "  --no-interpro: this flag will skip installation of interproscan, which uses a significant space on disk"
	echo "  --no-debian:   this flag will skip installation of Debian packages"
	echo "  --no-brew:     this flag will skip installation of Brew and its packages"
	echo "  --no-docker:   this flag will skip installation of Docker and its packages"
}

if [ -z "$1" ] ; then
  echo "Error: you must at least specify the target installation dir for pantagruel and other software:"
  usage
  exit 1
else
  if [[ "$1" == '-h' || "$1" == '--help' ]] ; then
    usage
    exit 0
  else
    SOFTWARE=$(readlink -f "$1")
  fi
fi
shift

installinterpro='true'
installdebian='true'
installbrew='true'
installdocker='true'
BINS=${HOME}/bin

while [[ ! -z "${@}" ]] ; do
  case "$1" in 
    --no-interpro)
      installinterpro='false' ;;
    --no-debian)
      installdebian='false' ;;
    --no-brew)
      installbrew='false' ;;
    --no-docker)
      installdocker='false' ;;
    *)
      BINS=$(readlink -f "$1") ;;
  esac
  shift
done

PATH=${PATH}:${BINS}

echo "Installation of Pantagruel and dependencies: ..."
echo ""

mkdir -p ${SOFTWARE}/ ${BINS}/

cd ${SOFTWARE}/

# Install Pantagruel pipeline and specific phylogenetic modules 
echo "get/update git repositories for Pantagruel pipeline"
if [ -d ${SOFTWARE}/pantagruel ] ; then
  cd ${SOFTWARE}/pantagruel && git pull && cd -
else
  git clone https://github.com/flass/pantagruel
fi
if [ $? != 0 ] ; then
  echo "ERROR: Unable to create/update git repository ${SOFTWARE}/pantagruel"
  exit 1
fi
cd ${SOFTWARE}/pantagruel && git submodule init && git submodule update && cd -
if [ $? != 0 ] ; then
  echo "ERROR: Unable to create/update git submodules"
  exit 1
fi
echo ""

#~ for repo in tree2 pantagruel ; do
  #~ if [ -d ${SOFTWARE}/${repo} ] ; then
    #~ cd ${repo} && git pull && cd -
  #~ else
    #~ git clone https://github.com/flass/${repo}
  #~ fi
  #~ if [ $? != 0 ] ; then
    #~ echo "ERROR: Unable to create/update git repository ${SOFTWARE}/${repo}"
    #~ exit 1
  #~ fi
#~ done
#~ echo ""

if [ "$installdebian" == 'true' ] ; then
  # basic dependencies, libs and standalone software, R and packages, Python and packages
  echo "get/update required Debian packages"
  #~ deppackages="git build-essential cmake gcc g++ lftp clustalo raxml libhmsbeagle1v5 mrbayes \
  #~ r-base-core r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt \
  #~ sqlite3 sqlite3-doc libmagick++-dev python python-scipy python-numpy python-biopython python-biopython-sql python-igraph cython \
  #~ mpi-default-bin mpi-default-dev mrbayes-mpi python-pip openjdk-8-jdk openjdk-8-jre cd-hit \
  #~ bioperl"
  # 'bioperl' should not be needed as only required for prokka and included in the Homebrew distribution; better not install it otherwise to avoid conflicts
  # 'cd-hit' no longer required in the pre-annotation step (prior to Prokka; uses mmseqs instead
  # 'libdw*' are required for Interproscan to work properly on a Ubuntu 18.04.1 LTS system - however that should ideally come with the Interproscan java executable
  # 'libgsl*' are required for proper functioningof MASH, even though it should rely on its brew dependency
  # 'openjdk-11*' are required for Interproscan since release 5.37-76.0; this package has to be obtained from PPA ppa:openjdk-r/ppa (otherwise the default repository will provide an java-11-openjdk packageg that installs JDK 10 [as of 26/09/2019])
  deppackages="git build-essential cmake gcc g++ gfortran lftp clustalo raxml libhmsbeagle1v5 mrbayes \
  r-base-core r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt r-cran-phytools r-cran-rcolorbrewer \
  sqlite3 sqlite3-doc libmagick++-dev python python-scipy python-numpy python-biopython python-biopython-sql python-igraph cython \
  mpi-default-bin mpi-default-dev mrbayes-mpi python-pip openjdk-11-jdk openjdk-11-jre \
  libdw1 libdw-dev libgsl23 libgsl-dev flex bison libgmp3-dev"
  if [ "$installbrew" == 'true' ] ; then
    deppackages="$deppackages linuxbrew-wrapper"
  else
    # GNU parallel should be installed via Homebrew as a dependency of Prokka;
	# without Homebrew, install GNU parallel via Debian
    deppackages="$deppackages parallel"
  fi
  if [ "$installdocker" == 'true' ] ; then deppackages="$deppackages docker\.io" ; fi
  sudo add-apt-repository -y ppa:openjdk-r/ppa
  sudo apt update
  sudo apt install -y $deppackages
  if [ $? != 0 ] ; then
    echo "ERROR: could not install all required Debian packages:"
    apt list $deppackages
    echo "ERROR: abort installation of Pantagruel"
    exit 1
  fi
  echo ""
fi

# install Python package BCBio.GFF
sudo -H pip install bcbio-gff
# install Python package bioscripts.convert
sudo -H pip install bioscripts.convert

# install R packages not supported by Debian
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep topGO | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
EOF
fi
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep phytools | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
install.packages('phytools', repos='${CRANmirror}')
EOF
fi
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep pvclust | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
install.packages('pvclust', repos='${CRANmirror}')
EOF
fi
echo "library('topGO')" | R --vanilla &> /dev/null
#~ checkexec "Could not install R package 'topGO'"
if [ $? != 0 ] ; then
 echo "WARNING: Could not install R package 'topGO'; testing of GO term enrichment in clade-specific gene sets will not be available"
fi
echo "library('pvclust')" | R --vanilla &> /dev/null
#~ checkexec "Could not install R package 'pvclust'"
if [ $? != 0 ] ; then
 echo "WARNING: Could not install R package 'pvclust'; confidence values on accessory genome-based clustering of genomes will not be computed"
fi
echo ""

if [ "$installbrew" == 'true' ] ; then
  # first try running brew - if never used before, this will lead to INTERACTIVE setup
  brew
  brewdir='/home/linuxbrew/.linuxbrew'
  if [ -d ${brewdir} ] ; then
    echo "Found that Homebrew recipes are installed at: ${brewdir}"
  else
    brewdir="${HOME}/.linuxbrew"
    if [ -d ${brewdir} ] ; then
      echo "Found that Homebrew recipes are installed at: ${brewdir}"
    else
      echo "Could not find the folder where Homebrew recipes are installed (not in the usual location /home/linuxbrew/.linuxbrew/ or ${HOME}/.linuxbrew/)"
      echo "  This may be all right, as long as this folder exists and that the locations of its relevant subfolders are documented in the environment variables \${PATH}, \${MANPATH} and \${INFOPATH}."
      brewdir=''
    fi
  fi
  if [ ! -z "${brewdir}" ] ; then
    # Configure Linuxbrew
    if [ ! -e ${HOME}/.bash_profile ] ; then
      # because when .bash_profile exists, it superseeds .profile, which won't be loaded automatically
      # thus if creating .bash_profile, restore the loading of .profile
      echo ". \"\$HOME/.profile\"" > ${HOME}/.bash_profile
    fi
    if [[ -z "$(grep PATH ${HOME}/.bash_profile | grep '\.linuxbrew')" ]] ; then
      echo "export PATH=\"${brewdir}/bin:\$PATH\"" >> ${HOME}/.bash_profile
      editedprofile=true
    fi
    if [[ -z "$(grep MANPATH ${HOME}/.bash_profile | grep '\.linuxbrew')" ]] ; then
      echo "export MANPATH=\"${brewdir}/share/man:\$MANPATH\"" >> ${HOME}/.bash_profile
      editedprofile=true
    fi
    if [[ -z "$(grep INFOPATH ${HOME}/.bash_profile | grep '\.linuxbrew')" ]] ; then
      echo "export INFOPATH=\"${brewdir}/share/info:\$INFOPATH\"" >> ${HOME}/.bash_profile
      editedprofile=true
    fi
    source ${HOME}/.bash_profile
  fi
  if [[ -z "$(echo ${PATH} | grep '\.linuxbrew')" || -z "$(echo ${MANPATH} | grep '\.linuxbrew')" || -z "$(echo ${INFOPATH} | grep '\.linuxbrew')" ]] ; then
    echo "Error: the environment variables \${PATH}, \${MANPATH} and \${INFOPATH} do not contain the path to the Homebrew installation"
    if [ -z "${brewdir}" ] ; then
      echo "definition of these variables could not be modified (by editing and sourcing the file ~/.bash_profile because the location of the Homebrew installation is not known"
    fi
    echo "exit now"
    exit 1
  fi 
  
  # install Prokka using brew
  if [[ -z "$(brew list prokka 2> /dev/null)" ]] ; then
    brew update
    brew doctor
    brew install brewsci/bio/prokka
    checkexec "Could not install Prokka using Brew" "Succesfully installed Prokka using Brew"
  else
    echo "found Prokka already installed with Brew:"
    brew info prokka
  fi
  echo ""
  
  # install MMSeqs using brew
  if [[ -z "$(brew list mmseqs2 2> /dev/null)" ]] ; then
    brew doctor
    brew install mmseqs2
    checkexec "Could not install MMSeqs using Brew" "Succesfully installed MMSeqs using Brew"
  else
    echo "found mmseqs2 already installed with Brew:"
    brew info mmseqs2
  fi  
  
  # install MASH using brew
  if [[ -z "$(brew list mash 2> /dev/null)" ]] ; then
    brew doctor
    brew install brewsci/bio/mash
    checkexec "Could not install MASH using Brew" "Succesfully installed MASH using Brew"
  else
    echo "found mash already installed with Brew:"
    brew info brewsci/bio/mash
  fi
  echo ""
fi

# fetch Pal2Nal script
if [ ! -x ${BINS}/pal2nal.pl ] ; then
  if [ ! -x ${SOFTWARE}/pal2nal.v14/pal2nal.pl ] ; then
    echo "get pal2nal"
    wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
    tar -xzf ${SOFTWARE}/pal2nal.v14.tar.gz
    chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
    checkexec "Could not install pal2nal.pl" "Succesfully installed pal2nal.pl"
  fi
  ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
  checkexec "Could not link pal2nal.pl executable to ${BINS}/" "Succesfully linked pal2nal.pl executable to ${BINS}/"
else
  echo "found pal2nal.pl executable: ${BINS}/pal2nal.pl"
  if [ -L ${BINS}/pal2nal.pl ] ; then
    echo " linking to $(readlink -f ${BINS}/pal2nal.pl)"
  fi
fi
echo ""

# fetch MAD program
if [ ! -x ${BINS}/mad ] ; then
  if [ ! -x ${SOFTWARE}/mad/mad ] ; then
    echo "get MAD"
    wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
    unzip ${SOFTWARE}/mad2-2.zip
    chmod +x ${SOFTWARE}/mad/mad
    checkexec "Could not install MAD" "Succesfully installed MAD"
  fi
  ln -s ${SOFTWARE}/mad/mad ${BINS}/
  checkexec "Could not link MAD executable to ${BINS}/" "Succesfully linked MAD executable to ${BINS}/"
else
  echo "found MAD executable: ${BINS}/mad"
  if [ -L ${BINS}/mad ] ; then
    echo " linking to $(readlink -f ${BINS}/mad)"
  fi
fi
echo ""

# fetch LSD program
if [ ! -x ${BINS}/lsd ] ; then
  if [ ! -x ${SOFTWARE}/lsd_unix ] ; then
    echo "get LSD"
    wget https://github.com/tothuhien/lsd-0.3beta/releases/download/v0.3.3/lsd_unix
    chmod +x ${SOFTWARE}/lsd_unix
    checkexec "Could not install LSD" "Succesfully installed LSD"
  checkexec "Could not link LSD executable to ${BINS}/" "Succesfully linked LSD executable to ${BINS}/"
  fi
  ln -s ${SOFTWARE}/lsd_unix ${BINS}/lsd
else
  echo "found LSD executable: ${BINS}/lsd"
  if [ -L ${BINS}/lsd ] ; then
    echo " linking to $(readlink -f ${BINS}/lsd)"
  fi
fi
echo ""

if [ "$installdocker" == 'true' ] ; then
  # set up Docker group (with root-equivalent permissions) and add main user to it
  # !!! makes the system less secure: OK within a dedicated virtual machine but to avoid on a server or desktop (or VM with other use)
  if [[ -z "$(grep docker /etc/group)" ]] ; then
    sudo groupadd docker
  fi
  if [[ -z "$(grep docker /etc/group | grep ${USER})" ]] ; then
    sudo usermod -aG docker $USER
    newgrp docker << EOF
EOF
  fi
  checkexec "Could not set group 'docker' or let user '$USER' join it"
  
  
  # install ALE using Docker --- OK within a virtual machine
  if [[ -z "$(sudo docker images | grep boussau/alesuite)" ]] ; then
    sudo docker pull boussau/alesuite
    checkexec "Could not install ALE suite using Docker"
  else
    echo "found ALE suite already installed with Docker:"
    sudo docker images | grep boussau/alesuite
  fi
fi

if [ "$installinterpro" == 'true' ] ; then
  # install Interproscan
  iphost="ftp://ftp.ebi.ac.uk"
  iploc="pub/software/unix/iprscan/5/"
  
  currIPversion=$(interproscan --version 2> /dev/null | head -n 1 | sed -e 's/InterProScan version //')
  lastIPversion=$(lftp -c "open ${iphost} ; ls -tr ${iploc}/[0-9].[0-9][0-9]*-[0-9][0-9]*.[0-9] ; quit" | tail -n 1 | awk '{print $NF}')
  if [[ -z "${currIPversion}" || "${currIPversion}" != "${lastIPversion}" ]] ; then
   echo "get Interproscan ${lastIPversion}"
   if [[ ! -x ${SOFTWARE}/interproscan-${lastIPversion}/interproscan.sh ]] ; then
     ipsourceftprep="${iphost}/${iploc}/${lastIPversion}/"
     ipsourcefile="interproscan-${lastIPversion}-64-bit.tar.gz"
     if [[ ! -e ${SOFTWARE}/${ipsourcefile} ]] ; then
       wget ${ipsourceftprep}/${ipsourcefile}
     fi
     if [[ ! -e ${SOFTWARE}/${ipsourcefile}.md5 ]] ; then
       wget ${ipsourceftprep}/${ipsourcefile}.md5
     fi
     dlok=$(md5sum -c ${ipsourcefile}.md5 | grep -o 'OK')
     ip=3
     # allow 3 attempts to downlaod the package, given its large size
     while [[ "$dlok" != 'OK' && $ip -gt 0 ]] ; do 
       ip=$(( $ip - 1 ))
       echo "dowload of ${ipsourcefile} failed; will retry ($ip times left)"
       wget ${ipsourceftprep}/${ipsourcefile}
       wget ${ipsourceftprep}/${ipsourcefile}.md5
       dlok=$(md5sum -c ${ipsourcefile}.md5 | grep -o 'OK')done
     done
     if [[ "$dlok" != 'OK' ]] ; then
       echo "ERROR: Could not dowload ${ipsourcefile}"
       exit 1
     fi
     tar -pxvzf ${SOFTWARE}/${ipsourcefile}
     checkexec "Could not uncompress Interproscan successfully"
   else
     echo "found InterProScan executable:"
     ls ${SOFTWARE}/interproscan-${lastIPversion}/interproscan.sh
   fi
   currIPversion=${lastIPversion}
   # make sure InterProscan uses the correct version of Java (1.8)
   wrongjava=$(${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh | grep -o 'Java version .* is required to run InterProScan.')
   if [ ! -z "${wrongjava}" ] ; then
     needjava=$(echo ${wrongjava} | sed -e 's/Java version \(.*\) is required to run InterProScan./\1/')
     currjava=$(readlink -f `which java`)
     goodjava=$(echo ${currjava} | sed -e "s/java-[0-9]\+-/java-${needjava}-/")
     echo "default version of Java is wrong: ${currjava}"
     echo "will find executable for desired version: ${needjava}"
     if [ ! -e ${goodjava} ] ; then
       goodjava=$(echo ${currjava} | sed -e "s/java-[0-9]\+-/java-${needjava}.0-/")
     fi
     if [ ! -e ${goodjava} ] ; then
       goodjava=$(ls -d $(dirname $(dirname $(dirname ${currjava})))/*${needjava}*/bin/java)
     fi
     if [ -z ${goodjava} ] ; then
       echo "ERROR: Could not find the required version of java for InterProScan:"
       ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh
       exit 1
     else
       echo "found version of java suitable for InterProScan: ${goodjava}"
     fi
      #~ export PATH=$(dirname ${goodjava}):$PATH
      #~ ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh
     #~ checkexec "Interproscan test was not successful"
     #~ if [[ -z "$(grep JAVA4INTERPROSCAN ${HOME}/.bashrc )" ]] ; then
       #~ echo "export JAVA4INTERPROSCAN=${goodjava}" >> ${HOME}/.bashrc
       #~ checkexec "Could not store JAVA4INTERPROSCAN in .bashrc"
       #~ editedrc=true
     #~ fi
     echo "will edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
     if [ ! -z "$(grep 'JAVA=\$(type -p java)' ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh)" ] ; then
       mv ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original && \
        sed -e "s#^JAVA=\$(type -p java)#JAVA=${goodjava}#g" ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original > ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
        chmod +x ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
        echo "fixed definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh (original script saved as '${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original')"
       checkexec "Could not edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
     else
       # the file was already modified, be more agressive with sed substitution
       mv ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 && \
        sed -e "s#^JAVA=.\+#JAVA=${goodjava}#g" ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 > ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
        chmod +x ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && rm ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 && \
        echo "fixed definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
       checkexec "Could not edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
     fi     
   fi
   ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh -i ${SOFTWARE}/interproscan-${currIPversion}/test_proteins.fasta -f tsv
   checkexec "Interproscan test was not successful"
   # link the exec file
   rm -f ${BINS}/interproscan
   ln -s ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${BINS}/interproscan
   checkexec "Failed to link Interproscan executable to ${BINS}/" "Succesfully linked Interproscan executable to ${BINS}/"
  else
   echo "found up-to-date version of Interproscan at $(ls -l `which interproscan` | awk '{print $NF}')"
  fi
fi

if [[ -z "$(grep 'export PATH=' ${HOME}/.bashrc | grep ${BINS})" ]] ; then
  echo "export PATH=\$PATH:${BINS}" >> ${HOME}/.bashrc
  editedrc=true
fi

rm -f ${BINS}/pantagruel
ln -s ${SOFTWARE}/pantagruel/pantagruel ${BINS}/
checkexec "Failed to link pantagruel executable to ${BINS}/" "Succesfully linked pantagruel executable to ${BINS}/"

echo "Installation of Pantagruel and dependencies: complete"


if [ "${editedprofile}" == 'true' ] ; then
  echo "the file '~/.bash_profile' was edited to set environment variables; please reload it with \`source ~/.bash_profile\` or reset your session for changes to take effect" 
fi
if [ "${editedrc}" == 'true' ] ; then
  echo "the file '~/.bashrc' was edited to set environment variables; please reload it with \`source ~/.bashrc\` or reset your session for changes to take effect" 
fi
