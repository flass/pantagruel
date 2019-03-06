This is a minimal dataset for the purpose of quick testing the pipeline.

Use the following commands:

```sh
# install
mkdir ./pantagruel_pipeline/
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
cd ..
./pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline

# run pipeline on a list of 10 genome accessions (9 RefSeq + 1 GenBank) to download from NCBI FTP
# + 1 genome to annotate (using the contigs from assembly GCF_900113725.1)
./pantagruel_pipeline/pantagruel/pantagruel -d testPTGdatabase -r ./ -f PANTAGFAM -I e.mail@institu.ti.on -L ./pantagruel_pipeline/pantagruel/data/NCBI_Assembly_accession_ids_test_10Brady -a ./pantagruel_pipeline/pantagruel/data/custom_genomes init
./pantagruel_pipeline/pantagruel/pantagruel -i ./testPTGdatabase/environ_pantagruel_testPTGdatabase.sh all
```
