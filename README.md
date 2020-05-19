![license] ![dockerhup-pulls]

# *Pantagruel*: a bioinformatic pipeline for the inference of gene evolution scenarios in bacterial pangenomes.

## News

### 23/03/2019: docker image on DockerHub! :whale:

The worries of installing Pantagruel are in the past! now you just have to download the docker iamge from the [Dockerhub] repository!

```sh
docker pull flass/pantagruel-dep:master-latest
```

See details [here](https://github.com/flass/pantagruel/blob/master/INSTALL.md#building-the-docker-image).
Note that the image is quite large (2.05 GB) so make sure you have the relevant space on your machine (and on the filesystem partition that hosts the docker client). Also, InterProScan (used in task 04) is NOT included in docker image, see [here](https://github.com/flass/pantagruel/blob/master/INSTALL.md#interproscantask-04-not-included-in-docker-image).

### 17/12/2019: :whale: Dockerfile
The [Dockerfile] (a 'recipe' to generate a Docker image) is now ready!  
Use it to generate a Docker image containing all Pantagruel dependencies; then use that image to run `pantagruel` commands within the environment provided by the container, as explained [here](https://github.com/flass/pantagruel/blob/master/INSTALL.md#the-container-worryless-way).  
NB: this may be used to generate a Singularity container as well, even though it was not tested yet.

### 26/03/2019: now on bioRxiv! 
The preprint describing the software and underlying methods is available on bioRxiv preprint server at:  
[https://www.biorxiv.org/content/10.1101/586495v3](https://www.biorxiv.org/content/10.1101/586495v2).

## Reference

If using this software, please cite:  
**Lassalle F, Veber P, Jauneikaite E, Didelot X. Automated Reconstruction of All Gene Histories in Large Bacterial Pangenome Datasets and Search for Co-Evolved Gene Modules with Pantagruel.” *bioRxiv* 586495.**  
**[doi: 10.1101/586495](doi.org/10.1101/586495)**


## Aim and description
*Pantagruel* provides an all-in-one software solution to reconstruct the complex evolutionary process of diversification of bacterial genomes.  

From a dataset of bacterial genomes, builds a database that describes the homology structure of all genes in the dataset -- the pangenome. With that, *Pantagruel* will first generate two key outputs:  
- a **reference (species) tree**, depicting the main signal of evolutionary relationships between **input genomes**;  
- **gene trees**, depicting the evolutionary relationships between gene sequences within a family of homologous genes.   

![pipeline1]

A **scenario of gene evolution** is then inferred for each gene family in the dataset by reconciling the topology of gene trees with the reference species tree in a probabilistic framework.  

Such scenario describes the likely events of gene **duplication, horizontal transfer and loss** (DTL model) that marked the gene family history. These events are annotated on the branch of the gene tree and the of the reference tree, and make their history consistent.
From these annotations, one can derive the **history of gain and loss** of this gene family over the reference tree of species, and follow the diversification of *gene lineages* within and across *genome lineages*.  

Gene tree/species tree reconciliation methods rely on the best available information to infer such scenarios, as they account for the phylogeny of genes; a probablilistic method is chosen to quantify the statistical support for inferences, in face of the large space of possible scenario and to account for the uncertainty in the input gene phylogeny.  

While probablilistic reconciliation methods are computationally costly, this pipeline uses innovative phylogenetic apporoaches based on the reduction of gene trees to their informative backbone, allowing their use in a resonable time on **datasets of 1,000+ bacterial genome** and covering **multiple species**.

![pipeline2]


Gene history data are then gathered in the database, which provides a way to:  
- quantify gene-to-gene **co-evolution**, i.e. a gene association score based on the evolutionary events (notably HGT) they shared (quantified at the gene lineage level rather than the whole gene family);  
- classify genes into **orthologous clusters** based on the gain/loss scenarios, from which one can define **clade-specific gene sets**.  

Two version of the pipeline are distributed:  

- a script version, which source code is adaptable and can be deployed on high-performance computing (HPC) "cluster" Linux systems;  

- (in development) a pre-compiled Docker image that can be deployed on pretty much any platform, including swarms of virtual machines (VMs). Future versions will be implemented using Philippe Veber's [Bistro](https://github.com/pveber/bistro) framework.

See below for instruction on software [installation] and [usage]. Impatients can go directly to running [examples].

--------------------

## Using Pantagruel

### Tasks

#### Initiating a Pantagruel database

The pipeline can be run using a single interface to deploy the several arms of the pipeline.  
It first requires to initiate the *Pantagruel* database, i.e. giving it a name, creating the base file structure, defining main options.
The generic syntax at the init stage is as follows, passing the key parameters using options:  
```sh
pantagruel -d db_name -r root_dir {-A refseq_assembly_folder | -a custom_assembly_folder} [other options] init
```  
Alternatively, the various parameters can be set directly in a Pantagruel configuration file specified with option `-i`:
```sh
pantagruel -d db_name -r root_dir -i config_file init
```  
The configuration file `config_file` can be generated by editing **a copy** of the [template environment script]. Note than it is only safe to edit the top parameters.

#### Running the pipeline / Building a Pantagruel database

Then, the pipeline can be run step-by-step by performing a specified task:
```sh
pantagruel -i config_file TASK
```
with `TASK` to be picked among the following (equivalent digit/number/keywords are separated by a `|`):
```
  0|00|fetch|fetch_data
       fetch public genome data from NCBI sequence databases and annotate private genomes;
	   also run a quick genome-to-genome distance estimation with MASH
  1|01|homologous|homologous_seq_families
       classify protein sequences into homologous families
  2|02|align|align_homologous_seq
       align homologous protein sequences and translate alignemnts into coding sequences
  3|03|sqldb|create_sqlite_db
       initiate SQL database and load genomic object relationships
  4|04|functional|functional_annotations
       use InterProScan to functionally annotate proteins in the database, including with Gene Ontology and metabolic pathway terms
  5|05|core|core_genome_ref_tree
       select core-genome markers and compute reference tree
  6|06|genetrees|gene_trees
       compute gene tree
  7|07|reconciliations
       compute species tree/gene tree reconciliations
  8|08|specific|clade_specific_genes
       classify genes into orthologous groups (OGs) and search clade-specific OGs
  9|09|coevolution
       quantify gene co-evolution and build gene association network

```  

Alternatively, several tasks ca be run at once by providing a space-separated string of tasks identifiers:  
```sh
pantagruel -i config_file TASK1 TASK2 ...
```
Finally, it is possible to run the whole pipeline at once, simply perform the `all` task:
```sh
pantagruel -i config_file all
```  
#### Dependencies between subsequent tasks 

Note there are **dependencies between tasks**, which must be carried on mostly sequentially:  
- *00, 01, 02, 03 tasks* each strictly depend on the previous step: 00 -> 01 -> 02 -> 03  
- *functional annotation task 04* is optional - though highly recomended - and depends on the previous task 01: 01 -> 04 
- *reference tree task 05* depends on previous *task 03* (and thus the previous ones)  
- *gene trees task 06* only depends on the previous *task 03* (and thus the previous ones) **IF** the `-c|--collapse` option is **NOT** used: 03 -> 06  
- however, if the `-c` option is specified, *task 06* (specifically step 6.4 when in HPC mode) is also dependent on *task 05*: 03 + 05 -> 06  
- *gene tree/species tree reconciliation task 07* strictly depends on the previous steps: 05 + 06 -> 07  
- *orthologous group clustering task 08* depends on previous *reconciliation step 07*: 07 -> 08  
- *co-evolution network task 09* depends on previous *reconciliation task 07*: 07 -> 09  
- but if run after *task 08*, an additional version of the co-evolution network will be made by collapsing the full network, grouping gene nodes by orthologous group: 07 + 08 -> 09

So all in all, you're better off running all the tasks sequentially, for instance using `pantagruel all`.

Importantly, it is recomended to use the alternative interface provided by HPC scripts to run intensive tasks `06` and `07` on high-performance computer (HPC) clusters for efficient and time-tractable computing; see [HPC scripts doc] below.

### Usage example  
Here is a standard examples of using `pantagruel` program.

First, to create a *new* database, we need to run the `init` task. To pass the key parameters, including where to create the database and its name, we will be using options:
```sh
pantagruel -d databasename -r /root/folder/for/database -f PANTAGFAM -I e.mail@institu.ti.on \
 -A /folder/of/public/genome/in/RefSeq/format init
```  
or

```sh
pantagruel -d databasename -r /root/folder/for/database -f PANTAGFAM -I e.mail@institu.ti.on \
 -L list_of_RefSeq_accession_ids -a /folder/of/custom/genomes init
```  
Then, to actually run the pipeline, we will execute the subsequent tasks.  
At this stage, no options need to (or can) be specified trough the command line, as all parameters are already defined
following the database intitiation stage (see above) and were stored in a configuration file. You will now simply have to specify where to find this configuration file with the `-i` option.
Unless you moved it, the configuration file should be where it has been created automatically,
at `${root_dir}/${db_name}/environ_pantagruel_${db_name}.sh`, with `${db_name}` and `${root_dir}` the arguments of `-d` and `-r` options on the `pantagruel init` call.  
So in our case, to execute the first three tasks, up to gene family sequence alignement, you can type the following command:  
```sh
pantagruel -i /root/folder/for/database/databasename/environ_pantagruel_databasename.sh fetch homologous align
```  
Note that this config file can be edited in-between tasks, for instance to change the location of key input files that you moved, or to tweak paramters - however this may cause issues in task dependencies (see above).

Please note that Pantagruel is still under active development and may evolve rapidly to fix bugs and solve issues. 
It is thus strongly recomended to update this software regularly using `git pull && git submodule update` in the `pantagruel/` git repository folder.
- TIP!!: once the database environment is loaded (by doing `source environ_pantagruel_databasename.sh`), you can use the alias `panup` for easy software updates.
     
**If Pantagruel is updated** in-between the running of tasks of a database project, 
it is higly recommended to **run the `init` task again (with the same options)** before proceeding to the next steps, 
to ensure the environment defined by the configuration file is compatible with the software. 
Regenerating the configuration file *will not erase other data*.

A simple way to regenerate the configuration file under the same parameters is to use the command:"
```sh
pantagruel -i previous_config_file --refresh init
```

### Input data

Note that for the sake of computing evolutionary analyses that have any meaning at all, Pantagruel requires that you provide a minimum of **four** genomes in input - 
ideally much more, as Pantagruel can easily deal with several hundreds of genome. 
This minimum number 4 is to be split at the user's discretion between RefSeq-type and 'custom assemblies, through `-A`/`-L` and `-a` options, respectively.  

Here is a view of what input data passed to Pantagruel should look like.  

When specifying accession ids to be downloaded from NCBI Assembly FTP using `-L|--refseq_list` or `--refseq_list4annot`, the list file which path is given as the option's argument should look like this:  
```sh
# for instance if used the command `pantagruel -d databasename -r /root/folder -L ./list_of_RefSeq_accession_ids init`
cat list_of_RefSeq_accession_ids
# GCF_000156855
# GCF_001026115.1
# GCF_001088845.2
# ...
# as many rows as there are genomes to study
# Note the trailing '.z' (z = 1,2,3,...) indicating the assembly version is optional 
# and if it is provided, it WILL BE IGNORED. 
# This is beacause the LAST version of the accession will always be returned.
# If you really want to work with an outdated version, please download it yourself
# and feed it to Pantagruel using '-A' option.
```

When using assemblies downloaded from NCBI RefSeq using options `-A|--refseq_ass` or `--refseq_ass4annot`, the folder which path is given as the option's argument should have a content looking like this:  
```sh
# for instance if used the command `pantagruel -d databasename -r /root/folder -A ./assemblies init`

ls -AF assemblies/
# GCF_000156855.2_ASM15685v2/
# GCF_001026115.1_ASM102611v1/
# GCF_001088845.1_8727_5_60/
# ... 
# as many separate assembly folders as there are genomes to study

ls -AF assemblies/GCF_001026115.1_ASM102611v1/
# GCF_001026115.1_ASM102611v1_assembly_report.txt      GCF_001026115.1_ASM102611v1_protein.faa.gz
# GCF_001026115.1_ASM102611v1_assembly_stats.txt       GCF_001026115.1_ASM102611v1_protein.gpff.gz
# GCF_001026115.1_ASM102611v1_cds_from_genomic.fna.gz  GCF_001026115.1_ASM102611v1_rna_from_genomic.fna.gz
# GCF_001026115.1_ASM102611v1_feature_count.txt.gz     GCF_001026115.1_ASM102611v1_translated_cds.faa.gz
# GCF_001026115.1_ASM102611v1_feature_table.txt.gz     GCF_001026115.1_ASM102611v1_wgsmaster.gbff
# GCF_001026115.1_ASM102611v1_genomic.fna.gz           GCF_001026115.1_ASM102611v1_wgsmaster.gbff.gz
# GCF_001026115.1_ASM102611v1_genomic.gbff             annotation_hashes.txt
# GCF_001026115.1_ASM102611v1_genomic.gbff.gz          assembly_status.txt
# GCF_001026115.1_ASM102611v1_genomic.gff.gz           md5checksums.txt
# each folder contains a set of files with all the assembly and annotation information
```

When providing your on 'custom' genomes, the folder which path is given as the argument of options `a|--custom_ass`should  have a content looking like this:  

```sh
# for instance if used the command `pantagruel -d databasename -r /root/folder -a ./user_genomes init`

ls -AF user_genomes/
# contigs/
# annotation/
# strain_infos_databasename.txt

ls -AF user_genomes/contigs/
# B03_1.fasta
# D03_1.fasta
# F03_1.fasta
# ... 
```
There should be as many separate genomic FASTA files in `user_genomes/contigs/` as there are genomes to study.

```sh
ls -AF user_genomes/annotation/
# B03_1/
# C03_1/
# D03_1/
# ... 
```
There can be an arbitrary number of annotation folders, only those which names match a genomic FASTA files in contigs/ will be considered

```sh
ls -AF user_genomes/annotation/B03_1/
# Rhizobium_endolithicum_Q54.fna
# Rhizobium_endolithicum_Q54.ffn
# Rhizobium_endolithicum_Q54.faa
# Rhizobium_endolithicum_Q54.gff
```
The file names within each annotation folder do not matter, only their extensions.

```sh
cat user_genomes/strain_infos_databasename.txt
# assembly_id	genus	species	strain	taxid	locus_tag_prefix
# B03_1	Rhizobium	endolithicum	Q54	1335060	REQ54
# ...
```
There should be as many rows as FASTA files in the `contigs/` folder, plus one header row.
Fields should be **tab**-separated; the header must contain these field names: `assembly_id`, `genus`, `species`, `strain`, `taxid`, `locus_tag_prefix`; their order does not matter.
The values in `assembly_id` and `locus_tag_prefix` fields must be unique per row.
The `assembly_id` field value must match the begin of contig file name and match exactly the annotation folder name.
The characters `'-'` and `'_'` are forbiden in the `locus_tag_prefix` field.


### Options

Options are detailed here:  
```
# for Pantagruel task 0-9: 

  _only one mandatory option_

    -i|--initfile      path to Pantagruel configuration file.
                        this file is generated at init stage, from the specified options.
  _facultative runtime options_

    -F|--FORCE        (no value) FORCE mode: will erase any pre-existing main folder for the task
                        (default: off, pre-exisitance of a folder will result in an early error)

    -R|--resume       (no value) try and resume the task from previous run that was interupted
                        (available for tasks 04-07)

    -N|--threads      specify the number of threads to use for (some) parrallelizable tasks (defaults to the maximum available))
                        (available for tasks: 00-02, 04-08)
						Note 1: this does not apply to ALE/ecceTERA reconciliations,
						  which jobs always run sequentially due to their high memory needs.
						Note 2: this does not apply to tasks run through the HPC script,
						  which have their own interface to define number of used CPUs.

    -z|--compress     will try and compress result file on the go (especially bulky files that won't be used much later
                        This will induce possible extra decompressing/re-generating data computing time
                        when resuming a task run with -R; avoid using compression when likely to have to resume later
                        (for the moment only available for tasks: 00 i.e. 'fetch')

# for Pantagruel task init:

  _mandatory options_

    -d|--dbname        string. database name

    -r|--rootdir       path to root directory where to create the database; defaults to current folder
  
    It is also necessary to specify an input genome dataset! 
    This is possible via -a, -A or -L options, or a mixture of them.
  
  _facultative options_

    -i|--initfile      path to Pantagruel configuration file.
                         a file can be derived (i.e. manualy curated) from 'environment_pantagruel_template.sh' template.
                         Parameters values specified in this file will override other options.
                         Can also be combined alone with --refresh to update the software version used for an existing database.

    --refresh          (no value) Use in combination with the -i option above to simply refresh the configuration file.
                         (e.g. after an update of the software) the program will simply re-run the `pantagruel [options] init` command
                         that has been previously used to generate the config file; hence there is no need to repeat any other option.
                         (even -d and -r options can be omitted if `pantagruel -i config_file --refresh init` is used)
                         Note that when options had quoted string arguments, unpredictable behaviour might occur;
                         please verify the outcome in the regenerated config file.
						 New options can be added _after_ the --refresh option to change the value of environment variables in the config file:
						   to set non-default values if not already, e.g. turn on collapsing:
						     `pantagruel -i config_file --refresh -c init`
						   to revert to default values, e.g. turn off collapsing:
						     `pantagruel -i config_file --refresh -n init`
                           to do both, e.g. turn off collapsing and switch to use ecceTERA reconciliation method:
                             `pantagruel -i config_file --refresh -n -e ecceTERA init`
                         Currently, only short options (e.g., -e or -n, NOT --rec_method or --no_collapse) are supported after --refresh.
                         
 General configuration options:

    -I|--iam           database creator identity (e-mail address is preferred)

    -f|--famprefix     alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; defaults to 'PANTAG'
                         the chosen prefix will be appended with a 'P' for protein families and a 'C' for CDS families.
 
    --path_to_interproscan  path to the InterProScan executable script, or to the folder containing an executable file named \`interproscan\`,
                              which itself should link to the script \`interproscan.sh\` that is found in the InterProScan software archive.
                              Defaults to the empty string, meaning that Pantagruel will look for the \`interproscan\` command in the \$PATH.
                             Using this option is mandatory to execute task 04 when when calling \`pantagruel\` through the docker image.
							 
    -u|--update_from  the new pantagruel database will be an update from a template/reference pantagruel database found at this path.
                        Requirements:
                          1) The genome set covered by the template db should be included in the genome set of the new db
                          2) the gene family prefix tag provided by option -F should be the same
                        As a result of the update, all gene family identifiers will correspond to the same
                        families between the datasets (unless when they are unique to the new genomes).
 Input options:
 
    -T|--taxonomy      path to folder of taxonomy database flat files. Defaults to $rootdir/NCBI/Taxonomy_YYYY-MM-DD (suffix is today's date)
                         if this is not containing the expected file, triggers downloading the daily dump from NCBI Taxonomy at task 00

    -A|--refseq_ass    path to folder of source genome assembly __folders__ containing flat files formated like NCBI Assembly RefSeq (no default value).
                         The assembly folders (one per genome) can be obtained on https://www.ncbi.nlm.nih.gov/assembly,
                         by making a keyword search to select a set of assemblies and downloading results with options:
                           Source Database = 'RefSeq' and File type = 'All file types (including assembly-structure directory)'.
                         A simple archive 'genome_assemblies.tar' (as obtained from the NCBI website) can be placed in that folder.
                         If user genomes are also provided, these RefSeq assemblies will be used as reference for their annotation.

    --refseq_ass4annot folder path. Same principle as -A, but WILL NOT be used in the study, only as a reference to annotate user genomes
                         (defaults to combined value of -A and -L options)

    -L|--refseq_list   file path. Same as -A|--refseq_ass, but just provide a list of NCBI Assembly accession ids (file with one accession id per row)
                         Accession ids are of the form GCx_yyyyyyyyy.z with x = {A|F} for GenBank and RefSeq, respectively, and y an z are any digit.
                         These accessions will be fetched from the NCBI FTP site using lftp.
                         Note the LAST version of the accession will be returned, i.e. the trailing '.z' part of the accession id is ignored.
                         These assemblies saved into a folder named after the value of the option:
                           for instance, \`-L /path/to/assemblist\` will save assembly folders in /path/to/assemblist_assemblies_from_ftp/.

    --refseq_list4annot file path. Same principle as -L, but WILL NOT be used in the study, only as a reference to annotate user genomes
                        (defaults to combined value of -A and -L options)
 
    -a|--custom_ass  path to folder of user-provided genomes (no default value). The specified folder must contain:
                      _mandatory_ 
                       - a 'contigs/' folder, where are stored multi-FASTA files of genome assemblies (one file per genome,
                          with extension '.fa', '.fasta' or '.fas' ...). Fasta file names will be truncated by removing
                          the '.fa' string and everything occuring after) and will be retained as the assembly_id (beware 
                          of names redundant with RefSeq assemblies).
                       - a 'strain_infos_${databasename}.txt' TAB-delimited file describing the organism, with ${databasename} the value of option -d"
                          columns should be headed with these fields (replace quotes and semicolons by tabs!):"
                           'sequencing_project_id'; 'genus'; 'species'; 'strain'; 'taxid'; 'locus_tag_prefix'
                         'sequencing_project_id' must match the name of a contig file (e.g. 'seqProjID.fasta')
                         'locus_tag_prefix' must match the prefix of ids given to CDS, proteins and genome regions (contigs)
                         in potentially provided annotation files (see below).
                         Note that for ensuring compatibility with dependencies (namely BioPython Nexus alignment parser
                         and ALE reconciliation program), the characters '-' and '_' are forbidden in the 'locus_tag_prefix' field.
                      _optional_ 
                       - an 'annotation/' folder, where are stored annotation files: 
                         - one mandatory in GFF 3.0 file format (with a '.gff' extension);
                          and optionally, the following files (with consistent ids!!):
                         - one in GenBank flat file format (with a '.gbk' extension);
                         - one in Fasta format containing CDS sequences (with a '.ffn' extension).
                         - one in Fasta format containing matching protein sequences (with a '.faa' extension).
                         These four files are produced when using Prokka for annotation; if at least one of the .gbk, .ffn or .faa
                         are missing, all three will be derived from the .gff source. Each genome annotation file set must be stored
                         in a separate folder, which name must match a contig file (e.g. 'seqProjID/' for 'seqProjID.fasta').
                      NOTE: to ensure proper parsing, it is strongly advised that any provided annotation was generated with Prokka
                      NOTE: to ensure uniform annotation of the dataset, it is advised to let Pantagruel annotate the contigs (calling Prokka)
 
 --strain_info   path to an optional custom strain information file, provided in the same format as described above for 'strain_infos_\$\{databasename\}.txt'
                      This is only taken into account in the basence of option -a, i.e. only when using options -A or -L to specify input genomes
                      from GenBank/RefSeq (or genomes with compliant formats). This allows to override automated genome code generation.
					  
    -V|--env_var    quoted string of the form: 'variable1=value1[,variable2=value2[,...]]'.
                     Will add these variables to the configuration file so they can be exported to the environment during tasks.
                     Can be useful to define custom values of generic variables, e.g. "refgenus=Escherichia,seqcentre=Sanger_Institute"

 Output: core genome / reference phylogeny options:
  
    -s|--pseudocore  integer, float <=1.0 or string. The minimum number or fraction of genomes in which a gene family should be present
                       to be included in the pseudo-core genome, i.e. the gene set which alignments will be concatenated for reference tree search.
                       A non-numeric value will trigger an INTERACTIVE prompt for search of an optimal value at the begining of task 'core'.
                       Defaults to the total number of genomes (strict core genome set).

    -t|--reftree     Newick format tree file path. Specifies a reference tree for reconciliation and clade-specific gene analyses;
                       cancels the computation of tree from the concatenate of (pseudo-)core genome gene during task 'core'.

    --core_seqtype   {cds|prot} defines the type of sequence that will be used to compute the (pseudo-)core genome tree (default to 'cds')

    --pop_lg_thresh  real. Defines the threshold of branch length for delinating populations in the reference tree 
                       (default: 0.0005 for nucleotide alignment-based tree; 0.0002 for protein-based)

    --pop_bs_thresh  real. Defines the threshold of branch support for delinating populations in the reference tree (default: 80)
    
    --rooting        string. Defines the method to root the reference tree during task 5|core_genome_ref_tree. 
                       Possible values are 'treebalance', 'MAD' and 'outgroup:SPECIESCODELIST' (default: 'treebalance'),
                       - 'treebalance' uses the '-f I' algorthm of RAxML to root the tree towards an optimal balance of branch lengths
                          on either sides of the root;
                       - 'MAD' uses the minimal ancestor deviation method described in \"Tria, et al. (2017) Nat. Ecol. Evol. 1, 0193\".
                       - 'outgroup:SPECIESCODELIST' will root according tothe specified outgroup(s), with SPECIESCODELIST a comma-sperated list of species ids:
                            'outgroup:SPECIESCODE' for rooting with a single species
                            'outgroup:SPECIESCODE1,SPECIESCODE2,... for mutilple species (in which case their MRCA in  the tree will be the outgroup)
                          Species ids can be either valid genome assembly ids of the relevant input genomes (typically a NCBI Assembly accession id),
                          or internal genome identifiers that are specifically in the Pantagruel database but often match the relevant Uniprot organism code.
                          The mapping between genome accession ids and organism codes is given in the file '03.database/genome_codes.tab' generated during task 3.
                          To use codes, you may thus want run task 3 first, then run task init again with this option to regenerate the config file with 
                          the desired outgroup organism codes and only then run task 5.

    -S|--snp_aln        reduce the core-genome alignment to SNPs

 Output: gene trees / reconciliations options:
     
    -H|--submit_hpc  full address (hostname:/folder/location) of a folder on a remote high-performance computating (HPC) cluster server.
                       This indicate that computationally intensive tasks, including building the gene tree collection
                       ('genetrees') and reconciling gene tree with species tree ('reconciliations') will be run
                       on a HPC server (only Torque/PBS and LSF job submission systemd are supported so far).
                       [support for core genome tree building ('core') remains to be implemented].
                       Instead of running the computations, scripts for HPC cluster job array submission will be generated automatically.
                       Data and scripts will be transfered to the specified address (the database folder structure
                       will be duplicated there, but only relevant files will be synced). Note that job submission
                       scripts will need to be executed manually on the cluster server.
                       If set at init stage, this option will be maintained for all tasks. However, the remote address
                       can be updated when calling a specific task; string 'none' cancels the HPC behaviour.

    -c|--collapse    (no value) Enable collapsing the rake clades in the gene trees (strongly recomended in datasets of size >50 genomes).

    -n|--no_collapse    (no value) disable collapsing the rake clades in the gene trees"
                       (default; use this option in combination with -i --refresh to restore default behaviour when -c was used in previous runs)."

    -C|--collapse_par  quoted string. Specifies parameters for collapsing the rake clades in the gene trees.
                       A single-quoted, semicolon-delimited string containing variable definitions must be provided.
                       Default is equivalent to providing the following string:
                          'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'

    -e|--rec_method   {ALE|ecceTERA} choose the method to reconcile gene trees and the species tree.
                       ALE (default): a probabilistic method to sample gene Duplication, Transfer and Loss (DTL) scenarios
                        by amalgamating the likelihood of bayesian samples of trees (doi:10.1093/sysbio/syt003;doi:10.1093/sysbio/syt054).
                        The likelihood-based approach can be heavy in memory use (several 10GB for one gene family scenario) and computation time.
                        The option '-c' (to collapse gene trees prior to reconciliation) efficiently mitigates this issue as it generally reduces
                        the compute time to minutes is highly recommmended when the dataset size grows (>50 bacterial genomes).
                       ecceTERA: a parsimony method to sample gene DTL scenarios by amalgamating the likelihood of bayesian samples of trees
                        under a model and procedure similar to ALE (doi:10.1093/bioinformatics/btw105).
                        The parsimony apporach allows the use of this methods on large-scale datasets within a reasonable time and using little memory
                        without having to resort to gene tree collapsing with option '-c' (but using it is possible and would make reconciliation even faster).

    -g|--genefam_list Path to gene family list file. Resticts the computation of gene trees and all subsequent analyses to a list of gene families.
                        This impacts all task from 06 and forward. The list has to be one gene family identifier per line.
                        Gene family ids have to refer to existing ones in the database, and therefore can only be defined after the running of task 02.
                        It is therefore advised to first run the pipeline up to task 02 (or equally up to 05) without this option,
                        and then to to set this paramter for the downstream computations.
                        This can be done by editing the value of 'genefamlist' variable in the configuration file or by using:
                          pantagruel -i configfile --refresh -g genelist init  (note it is important that -g option be placed after the --refresh option)
                        Reverting to the exhaustive computation behavior can be done similarly by setting 'genefamlist' variable to an empty value or by using:
                          pantagruel -i configfile --refresh -g '' init
   
 Output: Gene co-evolution options:"

    -q|--max_event_age Older relative age on the species tree (real value between 0.0 = tips and 1.0 = root) under which events will considered to compute co-evolution scores"
                       and to build the gene co-evolution network. Deeper branches of the species tree are often long and aglomerate long evolutionary periods into one time point."
                       As a result, gene histories involving old events mapped to these deep branches will be more likely to correlate in an unspecific way."
                       Default value is 0.5, meaning events older than half the height of the ultrametric species tree are not considered for co-evolution scoring.

# for any Pantagruel command calls:

    -h|--help          print this help message and exit.
```  

### HPC scripts: submission of intensive tasks to high-performance computer clusters

Tasks `06|genetrees` and `07|reconciliations` are computationally intensive due to the use of Bayesian algorithm, and due to the sheer number of homologous gene families for which a tree and a an evolution scenario need to be computed.  
Thankfully, most academic institutions will nowadays give you access to a HPC cluster, that provides: 
 - high-efficience compute nodes for demanding tasks (reconciliation can be very memory hungry, up to above 100GB for complex gene family scenarios and large species trees);
 - an interface to submit many similar individual jobs as _arrays_ of jobs. The structure of data handled by Pantagruel - many gene families expecting the same computational treatment - lend themselves perfectly to this sort of computing infracstructure.  
 
It is therefore highly recomended to use HPC clusters to deal with these intensive tasks if you can and have a dense dataset (or want it done with quickly).

For this sake, Pantagruel package provides an alternative to the main interface, using shell scripts for submission of jobs to the HPC cluster. So far, only **Torque/PBS** and **IBM LSF** cluster systems are suported.
Tasks are broken down into steps, as every step within tasks need to be completed for all gene families (at least those you want to include in dowstream analyses) before things are carried forward.

Here is how to proceed:

1) first, you should refresh your configuration file so to include the HPC parameters through the `-H` option:  
```sh
pantagruel -i previous_config_file --refresh -H hpchost:/where/you/will/set/your/database init
```

2) then run task `06` as you would normally do:
```sh
pantagruel -i previous_config_file 06
```
This will generate the list of gene families for which trees are to be computed, then send the files to the designated location on the HPC server. 
This copy step will only work provided you can connect via `ssh` to that host and that it does not require specifying the remote user account or to type in a password, i.e. that the simple command `ssh hpchost` would get you logged in your account without prompting you for a password.
This is easily achieved i) using secure (e.g. RSA) key pairs to log-in between local and remote hosts (see [SSH documentation](https://help.ubuntu.com/community/SSH/OpenSSH/Keys)) and *using no passphrase* (which is more secure anyway) 
and ii) by using a SSH config file (located in our local home folder as `~/.ssh/config`, see [ssh_config documetation](http://manpages.ubuntu.com/manpages/trusty/man5/ssh_config.5.html#files)) to describe the log-in details (remote user name, etc.)

the task will then stop short and tell you what to do next:
```
please connect to remote host hpchost and execute the following scripts in order 
(waiting for completion of all array jobs submitted by one script before executing the next):
- pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh [OPTIONS] ptg_config_file
- pantagruel_pipeline_06.2_HPC_collapse_gene_trees.sh [OPTIONS] ptg_config_file
- pantagruel_pipeline_06.3_HPC_bayesian_gene_trees.sh [OPTIONS] ptg_config_file
- pantagruel_pipeline_06.4_HPC_replace_spe_by_pops.sh [OPTIONS] ptg_config_file
- pantagruel_pipeline_06.5_HPC_populate_db_collapsed_clades.sh [OPTIONS] ptg_config_file
then copy back ouput files and updated database file by syncing the root folder from remote host to this host
```
Of course, Pantagruel should be installed on the HPC host! 

Logging onto the HPC host, you should visit the folder where you copied your database, and open the database's configuration file to edit the value of the `$ptgrepo` environment variable so to reflect where the pantagruel git repository has been cloned on that host.  
Then, you can run the scripts sequentially (waiting for full completion in-between each step!) as indicated above.  
Note that parameters relevant to the HPC submission can be specified for these scripts using options (each script comes with its own set of default values, for instance regarding requested resources like max memory allowance on the compute node):  
```sh
ptgrepo=/where/you/cloned/pantagruel
ptgscripts=${ptgrepo}/scripts
# to see the options:
${ptgscripts}/pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh --help
# example of options, to specify that you will use the LSF system and request 32GB and 4 CPUs on each compute node, and 24h of maximum walltime use of the node:
${ptgscripts}/pantagruel_pipeline_06.1_HPC_full_ML_gene_trees.sh --mem 32 --ncpus 4 --wth 24 --hpctype 'LSF' ptg_config_file
```


-------------

## Installing Pantagruel and its dependencies

This bioinformatic pipeline relies on a quite a few other pieces of software.  
To install them, please follow the indications in the [INSTALL] page.  

Main options are to use either: 
- the automated [install_dependencies.sh] shell script (for Debian systems; tested on Ubuntu 18.04)
- (NEW!) the [Dockerfile] to generate a Docker image containing all dependencies; use that image to run `pantagruel` commmands within the environment provided by the container

-------------

![repas]

[license]:https://img.shields.io/github/license/flass/pantagruel
[dockerhup-pulls]:https://img.shields.io/docker/pulls/flass/pantagruel-dep.svg
[repas]: https://github.com/flass/pantagruel/blob/master/pics/Pantagruels_childhood.jpg
[pipeline1]: https://github.com/flass/pantagruel/blob/master/pics/extract_cluster_concat_spetree_MLgenetrees.png
[pipeline2]: https://github.com/flass/pantagruel/blob/master/pics/collapse_samplebackbones_reconcile_compare.png
[installation]: https://github.com/flass/pantagruel#installing-pantagruel-and-its-dependencies
[INSTALL]: https://github.com/flass/pantagruel/blob/master/INSTALL.md
[usage]: https://github.com/flass/pantagruel#using-pantagruel
[examples]: https://github.com/flass/pantagruel#usage-example
[template environment script]: https://github.com/flass/pantagruel/blob/master/scripts/pipeline/environ_pantagruel_template.sh
[install_dependencies.sh]: https://github.com/flass/pantagruel/blob/master/install_dependencies.sh
[Dockerfile]:  https://github.com/flass/pantagruel/blob/master/Dockerfile
[HPC scripts doc]: https://github.com/flass/pantagruel#hpc-scripts-submission-of-intensive-tasks-to-high-performance-computer-clusters
[Dockerhub]: https://hub.docker.com/repository/docker/flass/pantagruel-dep