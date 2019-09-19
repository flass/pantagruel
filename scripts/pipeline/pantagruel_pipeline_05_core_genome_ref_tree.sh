#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

checkptgversion
checkfoldersafe ${coregenome}

if [ "${resumetask}" == "true" ] ; then
  echo "will try and resume computation where it was last stopped; may skip/resume computing: core genome concatenated alignment, core ML tree search, core tree bootstrapping, core tree rooting with ${rootingmethod}"
fi

case "${coreseqtype}" in
  prot)
      alifastacodedir=${protalifastacodedir}
      raxmloptions="-n ${treename} -m PROTCATLGX -j -p 1753 -w ${coretree}"
      raxmloptionsG="-n ${treename} -m PROTGAMMALGX -j -p 1753 -w ${coretree}"
      if [[ -z "${poplgthresh}" || "${poplgthresh}" == 'default' ]] ; then
        poplgthresh=0.0002
      fi
      ;;
  cds)
      alifastacodedir=${cdsalifastacodedir}
      raxmloptions="-n ${treename} -m GTRCATX -j -p 1753 -w ${coretree}"
      raxmloptionsG="-n ${treename} -m GTRGAMMAX -j -p 1753 -w ${coretree}"
      if [[ -z "${poplgthresh}" || "${poplgthresh}" == 'default' ]] ; then
        poplgthresh=0.0005
      fi
      ;;
  *)
      echo "ERROR: incorrect core sequence type was specified: ${coreseqtype} (must be 'prot' or 'cds'); exit now"
      exit 1 ;;
esac
popstemconds="[('lg', '>=', ${poplgthresh}, ${poplgleafmul}), ('bs', '>=', ${popbsthresh})]"
withinpopconds="[('max', 'lg', '<=', ${poplgthresh}, -1)]"


###########################################
## 05. Core Genome Phylogeny (Species tree)
###########################################

## select psedo-core genome marker gene set
export coregenome=${ptgdb}/05.core_genome
mkdir -p ${coregenome}/

if [[ ! -z "${userreftree}" ]] ; then
  # use user-provided tree as reference tree
  if [[ ! -e "${userreftree}" ]] ; then
    echo "Error: cannot open user-provided tree '$userreftree' ; exit now."
    exit 1
  else
    # test if it is a valid Newick tree file
    R --vanilla --slave << EOF
    library('ape')
    nftree = '${userreftree}'
    tree = read.tree(nftree)
    if (is.null(tree)){
      cat(sprintf("cannot read '%s' as a Newick tree file\n", nftree))
      quit(status=3)
    }else{
      goodtipset = TRUE
      # test that tip labels are all correct
      nfgencode = '${database}/genome_codes.tab'
      genome.codes = as.character(read.table(nfgencode, sep='\t', header=F)[,2])
      if (length(genome.codes)!=length(tree[['tip.label']])){
        goodtipset = FALSE
        cat("wrong number of tips compared to genome code set\n")
      }
      extratips = setdiff(tree[['tip.label']], genome.codes)
      missdtips = setdiff(genome.codes, tree[['tip.label']])
      if (length(extratips)>0){
        goodtipset = FALSE
        cat("tips not in genome code set:\n")
        print(extratips, quote=F)
      }
      if (length(missdtips)>0){
        goodtipset = FALSE
        cat("tips missing with respect to genome code set:\n")
        print(missdtips, quote=F)
      }
      if (!goodtipset){
        cat(sprintf("incorrect tree tip set (according to '%s')\n", nfgencode))
        quit(status=3)
      }
      if (!is.null(tree[['node.label']])){
        if (is.rooted(tree)){
          quit(status=6) 
        }else{
          quit(status=5) 
      }}else{
        quit(status=4)
      }
    }
EOF
    treecheck=$?
    # evaluate output from tree check
    case ${treecheck} in
     3)
       echo "ERROR: ${userreftree} is not a valid tree file; exit now" 1>&2
       exit 1 ;;
     4)
       echo "'$userreftree' is a correct tree file but has no branch support"
       echo "will use its topology for reference tree; skip RAxML ML tree computation but will compute rapid bootstraps from (pseudo-)core-genome concatenate alignment"
       echo "emulate resuming RAxML computation pipeline after ML tree search (not that tree will be re-rooted according to set criterion: '$rootingmethod')"
       cp ${userreftree} ${nrbesttree}
       resumetask='true' ;;
     5)
       echo "'$userreftree' is a correct tree file, with branch supports but not rooted"
       echo "will use its topology and branch supports for reference tree; skip RAxML ML tree and bootsrap computations but will perform rooting according to set criterion: '$rootingmethod'"
       echo "emulate resuming RAxML computation pipeline after rapid bootstrap tree search"
       cp ${userreftree} ${nrbesttree}
       ln -s $(realpath --relative-to=$(dirname ${nrbiparts}) ${nrbesttree}) ${nrbiparts}
       resumetask='true' ;;
     6)
      echo "'$userreftree' is a correct rooted tree file with branch supports"
      echo "will use it as a reference tree; skip all RAxML computations from (pseudo-)core-genome concatenate alignment"
       cp ${userreftree} ${speciestree} ;;
     *)
      echo "something went wrong when evaluatiing input tree (tree check output value: ${treecheck}); exit now"
      exit 1 ;;
    esac
  fi
fi

echo "run non-interactively '$ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh' to record the gene family set."
if [[ -z "${pseudocoremingenomes}" ]] ; then	
     echo "Default: will use a strict core-genome gene set, i.e. genes present in a single copy in all the ${ngenomes} studied genomes."
	 ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh ${initfile} ${ngenomes}
else	
     echo "Will use a pseudo-core gene set, i.e. genes present in a single copy in at least the ${pseudocoremingenomes} out of the ${ngenomes} studied genomes."
	 ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh ${initfile} ${pseudocoremingenomes}
fi

if  [[ "${resumetask}" == "true" && -s ${speciestree} ]] ; then
  echo "found reference tree in file '${speciestree}' ; skip reference tree inferrence"
else
  ### compute reference tree from (pseudo-)core genome
  
  ## generate core genome alignment path list
  if [[ "${resumetask}" == "true" && -s ${pseudocorealn} ]] ; then
   echo "skip concatenating core genome alignment"
  else
   rm -f ${coregenome}/pseudo-coregenome_sets/${pseudocore}_${coreseqtype}_aln_list
   for fam in `cat ${coregenome}/pseudo-coregenome_sets/${pseudocore}_families.tab` ; do
    ls ${alifastacodedir}/$fam.codes.aln >> ${coregenome}/pseudo-coregenome_sets/${pseudocore}_${coreseqtype}_aln_list
   done
   # concatenate pseudo-core prot/CDS alignments
   python ${ptgscripts}/concat.py ${coregenome}/pseudo-coregenome_sets/${pseudocore}_${coreseqtype}_aln_list ${pseudocorealn} > ${ptglogs}/concat_core_genome.log
   checkexec "failed to produce concatenated (pseudo)core-genome alignment" "created concatenated (pseudo)core-genome alignment in file '${pseudocorealn}'"
   rm ${alifastacodedir}/*_all_sp_new
  fi

  ## compute species tree using RAxML
  mkdir -p ${ptglogs}/raxml ${coretree}/
  # define RAxML binary and options; uses options -T (threads) -j (checkpoints) 
  if [[ ! -z "$(grep -o avx2 /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-AVX2)" ]] ; then
    raxmlflav='-AVX2 -U'
  elif [[ ! -z "$(grep -o avx /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-AVX)" ]] ; then
    raxmlflav='-AVX -U'
  elif [[ ! -z "$(grep -o sse3 /proc/cpuinfo | head -n 1)" && ! -z "$(which raxmlHPC-PTHREADS-SSE3)" ]] ; then
    raxmlflav='-SSE3 -U'
  else
    raxmlflav=''
  fi
  raxmlbin="raxmlHPC-PTHREADS${raxmlflav} -T $(nproc)"

  # first check the alignment for duplicate sequences and write a reduced alignment with  will be excluded
  if [[ "${resumetask}" == "true" && -e ${pseudocorealn}.identical_sequences ]] ; then
   echo "skip identical sequence removal in core alignment"
  else
   echo "# call: $raxmlbin -s ${pseudocorealn} ${raxmloptions} -f c" > ${ptglogs}/raxml/${treename}.check.log
   $raxmlbin -s ${pseudocorealn} ${raxmloptions} -f c &>> ${ptglogs}/raxml/${treename}.check.log
   checkexec "failed to remove identical sequence in core alignment" 
   if [ -e "${pseudocorealn}.reduced" ] ; then
     echo "removed identical sequence in core alignment; reduced alignemnt stored in file '${pseudocorealn}.reduced'"
   else
     echo "no identical sequence was found in the core alignment"
   fi
   grep 'exactly identical$' ${coretree}/RAxML_info.${treename} | sed -e 's/IMPORTANT WARNING: Sequences \(.\+\) and \(.\+\) are exactly identical/\1\t\2/g' > ${pseudocorealn}.identical_sequences
   mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_identical_sequences.${treename}
  fi
  if [[ ! -z "${userreftree}" ]] ; then 
    # work on the full alignment as the user-provided tree is expected to bear all leaves
    coretreealn=${pseudocorealn}
  elif [ -e ${pseudocorealn}.reduced ] ; then
    # in priority work on the reduced alignment (if there is one)
    coretreealn=${pseudocorealn}.reduced
  else
    coretreealn=${pseudocorealn}
  fi

  # search ML tree topology on (reduced) alignment nuder CAT-based model and with -F option
  # this means no final thorough tree optimization is conducted under GAMMA-based model
  # and the output file is thus RAxML_result.* (not RAxML_bestTree.*) ; to avoid overwriting by next step, it is automatically renamed as RAxML_resultTopo.* 
  if [[ "${resumetask}" == "true" && -s ${nrbesttopo} ]] ; then
   # ML tree search already done
   echo "skip ML tree topology search"
  else
   ckps=($(ls -t ${coretree}/RAxML_checkpoint.${treename}.* ${coretree}/RAxML_result*.${treename} 2> /dev/null))
   if [[ "${resumetask}" == "true" && -z "${ckps}" ]] ; then
    if [[ "${ckps[0]}" == "${coretree}/RAxML_result.${treename}" ]] ; then
      echo "found best topology tree file '${ckps[0]}'"
    else
      # resume search from checkpoint
      echo "resume ML tree topology search from checkpoint ${ckps[0]##*.}"
      mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bestTree.${treename}.up_to_ckp${ckps[0]##*.}
      mv ${ptglogs}/raxml/${treename}.ML_topo.log ${ptglogs}/raxml/${treename}.ML_topo.log.up_to_ckp${ckps[0]##*.}
      echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -F -t ${ckps[0]}" > ${ptglogs}/raxml/${treename}.ML_topo.log
      $raxmlbin -s ${coretreealn} ${raxmloptions} -F -t ${ckps[0]} &>> ${ptglogs}/raxml/${treename}.ML_topo.log
    fi
   else
    # initial search
    echo "ML tree topology search"
    echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -F" > ${ptglogs}/raxml/${treename}.ML_topo.log
    $raxmlbin -s ${coretreealn} ${raxmloptions} -F &>> ${ptglogs}/raxml/${treename}.ML_topo.log
   fi
   mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_result.${treename}
   mv ${coretree}/RAxML_result.${treename} ${coretree}/RAxML_resultTopo.${treename}
   rm -f ${nrbesttopo}
   ln -s $(realpath --relative-to=$(dirname ${nrbesttopo}) ${coretree}/RAxML_resultTopo.${treename}) ${nrbesttopo}
   checkexec "failed search for ML tree topology" "ML tree topology search complete; best tree stored in file '${nrbesttopo}'"
  fi
  
  # now optimize the branch length and model parameters on ML tree topology under GAMMA-based model with -f e option
  # the output file is also RAxML_result.* ; for clarity and backward comaptibility reasons, it is automatically renamed as RAxML_bestTree.*
  if [[ "${resumetask}" == "true" && -s ${nrbesttree} ]] ; then
   # ML tree search already done
   echo "skip ML tree parameter & branch length search"
  else
    if [[ -s ${coretree}/RAxML_bestTree.${treename} ]] ; then
      echo "found best tree parameter & branch length file '${coretree}/RAxML_bestTree.${treename}'"
    else
      echo "ML tree parameter & branch length search under GAMMA-based model"
      echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptionsG} -f e -t ${nrbesttopo}" > ${ptglogs}/raxml/${treename}.ML_brlen.log
      $raxmlbin -s ${coretreealn} ${raxmloptionsG} -f e -t ${nrbesttopo} &>> ${ptglogs}/raxml/${treename}.ML_brlen.log
    fi
    mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bestTree.${treename}
    mv ${coretree}/RAxML_result.${treename} ${coretree}/RAxML_bestTree.${treename}
    rm -f ${nrbesttree}
    ln -s $(realpath --relative-to=$(dirname ${nrbesttree}) ${coretree}/RAxML_bestTree.${treename}) ${nrbesttree}
   checkexec "failed search for ML tree parameter & branch length" "ML tree parameter & branch length search complete; best tree stored in file '${nrbesttree}'"
  fi

  # compute ${ncorebootstrap} rapid bootstraps (you can set variable by editing ${envsourcescript})
  if [[ "${resumetask}" == "true" && -s ${nrbiparts} ]] ; then
   # tree with branch supports already provided bootstrap tree search already done
   echo "skip bootstraping and mapping to ML tree"
  else
   if [[ "${resumetask}" == "true" && -s ${coretree}/RAxML_bootstrap.${treename} ]] ; then
    nbs=$(wc -l ${coretree}/RAxML_bootstrap.${treename} | cut -d' ' -f1)
    if [ ${nbs} -ge ${ncorebootstrap} ] ; then
     # bootstrap tree search already done
     echo "skip bootstrapping"
    else
     # resume bootstrap tree search from previous run
     export ncorebootstrap=$(( ${ncorebootstrap} - ${nbs} ))
     echo "resume bootstraping after ${nbs} trees (${ncorebootstrap} left to compute)"
     mv ${coretree}/RAxML_bootstrap.${treename} ${coretree}/RAxML_bootstrap.${treename}.firstbs${nbs}
     mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bootstrap.${treename}.up_to_bs${nbs}
     mv ${ptglogs}/raxml/${treename}_bs.log ${ptglogs}/raxml/${treename}_bs.log.up_to_bs${nbs}
     echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -x 198237 -N ${ncorebootstrap}" > ${ptglogs}/raxml/${treename}_bs.log
     $raxmlbin -s ${coretreealn} ${raxmloptions} -x 198237 -N ${ncorebootstrap} &>> ${ptglogs}/raxml/${treename}_bs.log
    fi
   else
    # bootstrapping tree search from scratch
    echo "search bootstrap trees"
    echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -x 198237 -N ${ncorebootstrap}" > ${ptglogs}/raxml/${treename}_bs.log
    $raxmlbin -s ${coretreealn} ${raxmloptions} -x 198237 -N ${ncorebootstrap} &>> ${ptglogs}/raxml/${treename}_bs.log
   fi
   # mapping bootstraps to ML tree
   mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bootstrap.${treename}
   echo "map branch supports on ML tree"
   echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -f b -z ${coretree}/RAxML_bootstrap.${treename} -t ${nrbesttree}" > ${ptglogs}/raxml/${treename}_mapbs.log
   $raxmlbin -s ${coretreealn} ${raxmloptions} -f b -z ${coretree}/RAxML_bootstrap.${treename} -t ${nrbesttree} &>> ${ptglogs}/raxml/${treename}_mapbs.log
   mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_bipartitions.${treename}
   rm -f ${nrbiparts}
   ln -s $(realpath --relative-to=$(dirname ${nrbiparts}) ${coretree}/RAxML_bipartitions.${treename}) ${nrbiparts}
   checkexec "failed bootstrapping" "bootstrapping complete; bootstrap trees stored in file '${nrbiparts}'"
  fi

  # root ML tree
  if [[ "${resumetask}" == "true" && -s ${nrrootedtree} ]] ; then
   echo "skip tree rooting with '${rootingmethod}' method"
  else
   # variables automoatically defined in ${envsourcescript} :
   #~ export nrbesttree=${coretree}/RAxML_bestTree.${treename}
   #~ export nrbiparts=${nrbesttree/bestTree/bipartitions}
   #~ export nrrootedtree=${nrbiparts}.${rootingmethod/:/-}rooted
   if [[ "${rootingmethod}" == 'treebalance' ]] ; then
    # root ML tree with RAxML tree-balance algorithm
    nrbipartbranchlab=$(readlink -f $nrbiparts | sed -e 's/bipartitions/bipartitionsBranchLabels/')
    if [ ! -s ${nrbipartbranchlab} ] ; then
      echo "Error: could not detect tree file with brach support annotated on the branches ; exit now"
      exit 1
    fi
    echo "# call: $raxmlbin -s ${coretreealn} ${raxmloptions} -f I -t ${nrbipartbranchlab}" > ${ptglogs}/raxml/${treename}_root.log
    $raxmlbin -s ${coretreealn} ${raxmloptions} -f I -t ${nrbipartbranchlab} &>> ${ptglogs}/raxml/${treename}_root.log
    mv ${coretree}/RAxML_info.${treename} ${coretree}/RAxML_info_rootedTree.${treename}
    rm -f ${nrrootedtree}
    ln -s $(realpath --relative-to=$(dirname ${nrrootedtree}) ${coretree}/RAxML_rootedTree.${treename}) ${nrrootedtree}
   else
    if [[ "${rootingmethod}" == 'MAD' ]] ; then
     # root tree with MAD (Tria et al. Nat Ecol Evol (2017) Phylogenetic rooting using minimal ancestor deviation. s41559-017-0193 doi:10.1038/s41559-017-0193)
     R BATCH --vanilla --slave << EOF
     source('${ptgscripts}/mad_R/MAD.R')
     mad.rooting = MAD(readLines('${nrbiparts}'), 'full')
     write(mad.rooting[[1]], file='${nrbiparts}.${rootingmethod}rooted')
     save(mad.rooting, file='${nrbiparts}.${rootingmethod}rooting.RData}')
EOF
    rm -f ${nrrootedtree}
    ln -s $(realpath --relative-to=$(dirname ${nrrootedtree}) ${nrbiparts}.${rootingmethod}rooting) ${nrrootedtree}
    else
     if [[ "${rootingmethod:0:8}" == 'outgroup' ]] ; then
      # outgroup rooting; assume argument was of the shape 'outgroup:SPECIESCODE' or 'outgroup:SPECIESCODE1,SPECIESCODE2,SPECIESCODE3,...'
      outgroup=$(echo ${rootingmethod} | cut -d':' -f2)
      roottag=$(echo ${rootingmethod} | tr ':' '-' | tr ',' '-')_rooted
      
      python << EOF
import tree2
nftree = '${nrbiparts}'
tree = tree2.read_newick(nftree)
nfacc2code = '${database}/genome_codes.tab'
with open(nfacc2code, 'r') as facc2code:
  dacc2code = dict([tuple(line.rstrip('\n').split('\t')) for line in facc2code])

outg = '${outgroup}'.split(',')
# verify the presence of specified outgroup species in tree
knownspe = list(set(outg) & set(tree.get_leaf_labels()))
unknownspe = list(set(outg) - set(tree.get_leaf_labels()))
translatedspe = []
if unknownspe:
  print "found %d of the specified outgroup species are absent from tree '%s'"%(len(unknownspe), nftree)
  print "try tointerpret them as genome assembly ids; translate into appropriate genome/species code from '%s'"%nfacc2code
  for uks in unknownspe:
    if uks not in dacc2code:
      raise KeyError, "the outgroup id '%s' is also absent from genome assembly list"
    else:
      translatedspe.append(dacc2code[uks])

mrcaoutg = tree.mrca(knownspe+translatedspe)
tree.newOutgroup(mrcaoutg)
tree.write_newick('${nrbiparts}.${roottag}', comment='comment')
EOF
      #~ if [ -z $(grep ${outgroup} ${nrbesttree}) ] ; then
        #~ echo "Error, outgroup species '$outgroup' is absent from tree '${nrbesttree}'"
        #~ exit 1
      #~ fi
      #~ R --vanilla --slave << EOF
      #~ library('ape')
      #~ tree = read.tree('${nrbiparts}')
      #~ rtree = root(tree, outgroup='${outgroup}', resolve.root=T)
      #~ write.tree(rtree, file='${nrbiparts}.${roottag}')
#~ EOF
      rm -f ${nrrootedtree}
      ln -s $(realpath --relative-to=$(dirname ${nrrootedtree}) ${nrbiparts}.${roottag}) ${nrrootedtree}
     fi
    fi
   fi
  fi
  checkexec "failed rooting reference tree" "reference tree rooting complete; rooted tree stored in file '${nrrootedtree}'"

  # from here no skipping in resume mode
  if [ "${coretreealn}" == "${pseudocorealn}.reduced" ] ; then
    # reintroduce excluded species name into the tree
    python ${ptgscripts}/putidentseqbackintree.py --input.nr.tree=${nrbiparts} --ref.rooted.nr.tree=${nrrootedtree} \
     --list.identical.seq=${pseudocorealn}.identical_sequences --output.tree=${speciestree}
  else
    # only process tree format 
    python ${ptgscripts}/putidentseqbackintree.py --input.nr.tree=${nrbiparts} --ref.rooted.nr.tree=${nrrootedtree} \
     --output.tree=${speciestree}
  fi
  checkexec "failed re-introducing identical sequences into reference tree" "succesfully re-introduced identical sequences into reference tree; full reference tree stored in file '${speciestree}'"
fi

### make version of full reference tree with actual organism names
python ${ptgscripts}/code2orgnames_in_tree.py ${speciestree} $database/organism_codes.tab ${speciestree}.names
checkexec "failed name translation in reference tree" "reference tree name translation in complete; tree with organism names stored in file '${speciestree}.names'"

### generate ultrametric 'dated' species tree for more meaningfulgptghical representation of reconciliation AND to use the dated (original) version of ALE
## use LSD (v 0.3 beta; To et al. Syst. Biol. 2015) to generate an ultrametric tree from the rooted ML tree (only assumin different, uncorrelated clocks for each branch)
if [ -e "${pseudocorealn}.reduced" ] ; then
  alnlen=$(head -n1 ${pseudocorealn}.reduced | cut -d' ' -f2 )
elif [ -e "${pseudocorealn}" ] ; then
  #~ alnlen=$(head -n1 ${pseudocorealn} | cut -d' ' -f2 )
  alnlen=$(python << EOF
nfaln = '${pseudocorealn}'
alnlen = 0
nseq = 0
with open(nfaln, 'r') as faln:
  for line in faln:
    if line.startswith('>'):
      if nseq==0: nseq += 1
      else: break
    else:
      alnlen += (len(line) - 1)

print alnlen
EOF
)
else
  alnlen=500000 # arbitrary length similar to a 500 CDS concatenate
fi
lsd -i ${speciestree} -c -v 1 -s ${alnlen} -o ${speciestree}.lsd
checkexec "failed ultrametrization of reference tree" "ultrametrization of reference tree complete; ultrametric tree stored in file '${speciestree}.lsd'"

# export this ultrametric tree to Newick format (latest verison of LSD don't do it any more)
python << EOF
import sys
sys.setrecursionlimit(20000)
import tree2
nfin = '${speciestree}.lsd.nexus'
nfout = '${speciestree}.lsd.nwk'
t = tree2.read_nexus(nfin, returnDict=False, allLower=False, checkNewick=True)[0]
tree2.write_newick(t, nfout, comment=None)
EOF

### delineate populations of near-identical strains (based on tree with branch lengths in subs/site) and generate the population tree, i.e. the species tree withs population collapsed
python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -S ${speciestree} --threads=8 \
--pop_stem_conds="${popstemconds}" --within_pop_conds="${withinpopconds}"
checkexec "failed search of populations in reference tree" "search of populations in reference tree complete; population description stored in file '${speciestree%.*}_populations'"
nspepop=$(tail -n+3 ${speciestree%.*}_populations | wc -l)
echo "Found $nspepop disctinct populations in the species tree"


### extract sister clade pairs from the reference species tree for later clade-specific gene search
python << EOF
import tree2, sys
sys.setrecursionlimit(20000)
reftree = tree2.read_check_newick('${speciestree}')
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
checkexec "failed defining contrasting clades in reference tree" "defining contrasting clades in reference tree complete; clade description stored in file '${speciestree}_clade_defs'"

if [ -s ${indata}/all_assemblies_mash.dist ] ; then
    # if MASH was run during task 00
    Rscript ${ptgscripts}/plotmashdistcluster.r ${indata}/all_assemblies_mash.dist ${genomeinfo}/assembly_metadata/metadata.tab ${speciestree} ${database}/genome_codes.tab
fi
