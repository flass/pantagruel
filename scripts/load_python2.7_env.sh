#!/bin/bash

echo "\$PYTHONPATH=$PYTHONPATH"
python -c "import tree2"
while [ ${?} -gt 0 ] ; do
  echo "tree2 module was not detected."
  if [ ! -z "$1" ] ; then
    echo "try with path given through argument: $1"
	export PYTHONPATH=$1
    shift
  else
    echo "try to reload $HOME/.bashrc"
    source $HOME/.bashrc
	break
  fi
  python -c "import tree2"
done
echo ""

# use Anaconda to load required Python modules
# requires prior creation of the environment with:
# `conda create -n ptgenv python=2.7 pip ; conda activate ptgenv ; pip install scipy numpy biopython bcbio-gff Cython igraph psycopg2`
pyenvwarn="Warning: was not able to automatically load the Python environment ; the subsequent call to Python scripts may fail."
if [ ! -z "$(conda env list | grep ptgenv 2> /dev/null)" ] ; then
  conda activate ptgenv
elif  [ ! -z $(command -v module) ] ; then
  # use ad-hoc environment loading (deprecated: only known on Imperial College London HPC in its vintage form)
  #module load python
  module load anaconda3/personal
#  source activate env_python2
  source activate ptgenv
else
  echo ${pyenvwarn}
fi
if [ ${?} -gt 0 ] ; then
  echo ${pyenvwarn}
fi