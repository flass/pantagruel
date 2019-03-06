#!/bin/bash
echo "test installing pipeline"
mkdir ./pantagruel_pipeline/
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
cd ..
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline
echo "test running pipeline"
./pantagruel_pipeline/pantagruel/pantagruel -d testPTGdatabase -r ./ -f PANTAGFAM -I e.mail@institu.ti.on -L ./pantagruel_pipeline/pantagruel/data/NCBI_Assembly_accession_ids_test_10Brady -a ./pantagruel_pipeline/pantagruel/data/custom_genomes init
pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh all
