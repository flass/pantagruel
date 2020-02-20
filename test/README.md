## Testing `Pantagruel`

A minimal dataset is enclosed in the folder `data/` of the repository for the purpose of quick testing the pipeline.  
This dataset comprises a list of 10 genome accessions (9 from RefSeq + 1 from GenBank) to download from NCBI FTP + 1 'custom' genome (the bare contigs from assembly `GCF_900113725.1`) to test the annotate by `prokka` with a custom reference genus-specific database built from the other 10 genomes.

### Script-based install
To test the script-based install, please use the following commands:

```sh
# fetch pantagruel source code
mkdir ./pantagruel_pipeline/
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
cd ..
# install
./pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/

# run pipeline tests
pantagruel -d testPTGdatabase -r ./ -f PANTAGFAM -I e.mail@institu.ti.on -L ./pantagruel_pipeline/pantagruel/data/NCBI_Assembly_accession_ids_test_10Brady -a ./pantagruel_pipeline/pantagruel/data/custom_genomes init
pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh --refresh -c -g ./pantagruel_pipeline/pantagruel/data/test_genefam_list -q 0.75 init
pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh all
```

### Dockerfile-based install
To test the Dockerfile-based install, please use the following commands:

```sh
# fetch pantagruel source code
mkdir ./pantagruel_pipeline/
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
cd ..

# build the docker image for pantagruel dependencies
docker build -t panta etc
# install interproscan separately
pantagruel_pipeline/pantagruel/install_install_interproscan.sh pantagruel_pipeline/


# run pipeline tests 
docker run -u $UID -v $PWD:$PWD -w $PWD panta pantagruel -d testPTGdatabase -r ./ -f PANTAGFAM -I e.mail@institu.ti.on -L ./pantagruel_pipeline/pantagruel/data/NCBI_Assembly_accession_ids_test_10Brady -a ./pantagruel_pipeline/pantagruel/data/custom_genomes init
docker run -u $UID -v $PWD:$PWD -w $PWD panta pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh --refresh -c -g ./pantagruel_pipeline/pantagruel/data/test_genefam_list -q 0.75 init
docker run -u $UID -v $PWD:$PWD -w $PWD panta pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh all
```

Scripts are provided to run the testing part of the above in a more flexible way, see [test_install_run_pantagruel.sh](https://github.com/flass/pantagruel/blob/master/test/test_install_run_pantagruel.sh) and [test_docker_run_pantagruel.sh](https://github.com/flass/pantagruel/blob/master/test/test_docker_run_pantagruel.sh).
