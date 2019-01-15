#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}

if [ -z $3 ]
  pseudocoremingenomes=$3
fi

###########################################
## 05. Core Genome Phylogeny (Species tree)
###########################################

## select psedo-core genome marker gene set
export coregenome=${ptgdb}/05.core_genome
mkdir -p ${coregenome}/

if [[ -z ${pseudocoremingenomes} || "${pseudocoremingenomes}"=='REPLACEpseudocoremingenomes' ]] ; then
  echo "Error: 'pseudocoremingenomes' variable is not set; please run $ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh (interactive) to choose a sensible value."
  exit 1
fi 
## generate core genome alignment path list

rm -f ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list
for fam in `cat ${coregenome}/pseudo-coregenome_sets/${pseudocore}_families.tab` ; do
 ls ${protalifastacodedir}/$fam.codes.aln >> ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list
done
# concatenate pseudo-core prot alignments
python ${ptgscripts}/concat.py ${coregenome}/pseudo-coregenome_sets/${pseudocore}_prot_aln_list ${pseudocorealn}
rm ${protalifastacodedir}/*_all_sp_new

### compute species tree using RAxML
mkdir -p ${ptglogs}/raxml ${coretree}
# define RAxML binary and options; uses options -T (threads) -j (checkpoints) 
if [ ! -z $(grep -o avx2 /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-AVX2 -U'
elif [ ! -z $(grep -o avx /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-AVX -U'
elif [ ! -z $(grep -o sse3 /proc/cpuinfo | head -n 1) ] ; then
  raxmlflav='-SSE3 -U'
else
  raxmlflav=''
fi
raxmlbin="raxmlHPC-PTHREADS${raxmlflav} -T $(nproc)"
raxmloptions="-n ${treename} -m PROTCATLGX -j -p 1753 -w ${coretree}"

## first check the alignment for duplicate sequences and write a reduced alignment with  will be excluded
$raxmlbin -s ${pseudocorealn} ${raxmloptions} -f c &> ${ptglogs}/raxml/${treename}.check.
grep 'exactly identical$' ${coretree}/RAxML_info.${treename} | sed -e 's/IMPORTANT WARNING: Sequences \(.\+\) and \(.\+\) are exactly identical/\1\t\2/g' > ${pseudocorealn}.identical_sequences
rm ${coretree}/RAxML_info.${treename}
## first just single ML tree on reduced alignment
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} &> ${ptglogs}/raxml/${treename}.ML_run.log
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bestTree.${treename}

## RAxML, with rapid bootstrap
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -x 198237 -N 1000 &> ${ptglogs}/raxml/${treename}_bs.log
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bootstrap.${treename}
$raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -f b -z ${coretree}/RAxML_bootstrap.${treename} -t ${coretree}/RAxML_bestTree.${treename} 
mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bipartitions.${treename}

## root ML tree
# export nrbesttree=${coretree}/RAxML_bestTree.${treename}
# export nrbiparts=${nrbesttree/bestTree/bipartitions}
# export nrrootedtree=${nrbiparts}.${rootingmethod}rooted
if [[ "${rootingmethod}" == 'treebalance' ]] ; then
 # root ML tree with RAxML tree-balance algorithm
 $raxmlbin -s ${pseudocorealn}.reduced ${raxmloptions} -f I -t ${coretree}/RAxML_bipartitionsBranchLabels.${treename}  
 mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_rootedTree.${treename}
 ln -s ${coretree}/RAxML_rootedTree.${treename} ${nrrootedtree}
else
if [[ "${rootingmethod}" == 'MAD' ]] ; then
 # root tree with MAD (Tria et al. Nat Ecol Evol (2017) Phylogenetic rooting using minimal ancestor deviation. s41559-017-0193 doi:10.1038/s41559-017-0193)
 R BATCH --vanilla --slave << EOF
 source('${ptgscripts}/mad_R/MAD.R')
 mad.rooting = MAD(readLines('${nrbesttree}'), 'full')
 write(mad.rooting[[1]], file='${nrrootedtree}')
 save(mad.rooting, file='${nrbesttree}.${rootingmethod}rooting.RData}')
EOF
else
if [[ "${rootingmethod:0:8}" == 'outgroup' ]] ; then
 # outgroup rooting; assume argument was of the shape 'outgroup:SPECIESCODE'
 # verify the presence of outgroup species in tree
 outgroup=$(echo ${rootingmethod} | cut -d':' -f2)
 if [ -z $(grep ${outgroup} ${nrbesttree}) ] ; then
   echo "Error, outgroup species '$outgroup' is absent from tree '${nrbesttree}'"
   exit 1
 fi
 R BATCH --vanilla --slave << EOF
 library('ape')
 tree = read.tree('${nrbesttree}')
 root(tree, outgroup='${outgroup}')
EOF
 
fi ; fi ; fi
# reintroduce excluded species name into the tree + root as MAD-rooted species tree
export speciestree=${nrrootedtree}.full
python ${ptgscripts}/putidentseqbackintree.py --input.nr.tree=${nrbiparts} --ref.rooted.nr.tree=${nrrootedtree} \
 --list.identical.seq=${pseudocorealn}.identical_sequences --output.tree=${speciestree}

python ${ptgscripts}/code2orgnames_in_tree.py ${speciestree} $database/organism_codes.tab ${speciestree}.names

### generate ultrametric 'dated' species tree for more meaningfulgptghical representation of reconciliation AND to use the dated (original) version of ALE
## use LSD (v 0.3 beta; To et al. Syst. Biol. 2015) to generate an ultrametric tree from the rooted ML tree (only assumin different, uncorrelated clocks for each branch)
alnlen=$( head -n1 ${pseudocorealn}.reduced | cut -d' ' -f2 )
lsd -i ${speciestree} -c -v 1 -s $alnlen -o ${speciestree}.lsd

## delineate populations of near-identical strains (based on tree with branch lengths in subs/site) and generate the population tree, i.e. the species tree withs population collapsed
python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -S ${speciestree} --threads=8 \
--pop_stem_conds="[('lg', '>=', 0.0005), ('bs', '>=', 80)]" --within_pop_conds="[('max', 'lg', '<=', 0.0005, -1)]"
nspepop=$(tail -n+3 ${speciestree%.*}_populations | wc -l)
echo "Found $nspepop disctinct populations in the species tree"


## extract sister clade pairs from the reference species tree for later clade-specific gene search
python << EOF
import tree2
reftree = tree2.Node(file='${speciestree}')
nfout = "${speciestree}_clade_defs"
fout = open(nfout, 'w')
fout.write('\t'.join(['', 'clade', 'sisterclade'])+'\n')
k = 0
for node in reftree:
  if len(node.children) != 2: continue
  for child in node.children:
    if child.nb_leaves() <= 1: continue
    focchildlab = child.label()
    if not focchildlab:
      focchildlab = "clade%d"%k
      k += 1
    focchildleaflabset = ','.join(sorted(child.get_leaf_labels()))
    sischildleaflabset = ','.join(sorted(child.go_brother().get_leaf_labels()))
    fout.write('\t'.join([focchildlab, focchildleaflabset, sischildleaflabset])+'\n')

fout.close()
EOF
