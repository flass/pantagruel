#!/usr/bin/env Rscript
library(phytools)

combineNames = function(ass, metadata){
i = which(metadata[,1]==ass)
  if (length(i)==0){ 
    return(ass)
  }else{ if (length(i)>1){ 
    stop(sprintf("non-unique match to '%s' in metadata file", ass))
  }else{ 
    org = as.character(metadata[i,2])
  strain = as.character(metadata[i,4])
  if (nchar(strain)>0){
    nrorg = gsub(strain, '', org, ignore.case=T)
  }else{
    nrorg = org
  }
    ln = gsub(',', '-', gsub(' ', '_', paste(nrorg, strain)))
    return(paste(ass, ln, sep='_|_', collapse='_'))
  }}
}

plottreeheatmap = function(t, smat, pdfout, wid=45, hei=30){
  pdf(pdfout, width=wid, height=hei)
  ut = untangle(t, "read.tree")
  ruttl = rev(ut$tip.label)
  print(setdiff(ruttl, rownames(smat)))
  phylo.heatmap(t, (1.0 - smat[, ruttl]))
  dev.off()  
}

cargs = commandArgs(trailingOnly=T)
nfmashtriangle = cargs[1]
nfassmetadata = cargs[2]

if (length(cargs)>2){
  nftree = cargs[3]
  print(nftree)
  reftree = read.tree(nftree)
}else{
  reftree = NULL
}

if (length(cargs)>3){
  nftreelab2ass = cargs[4]
  print(nftreelab2ass)
  treelab2ass = read.table(nftreelab2ass, h=F, sep='\t', comment.char='')
}else{
  treelab2ass = NULL
}
print(nfmashtriangle)
mashdist = readLines(nfmashtriangle)

ntaxa = as.numeric(strsplit(mashdist[[1]], split='\t')[[1]][2])

ml = t(sapply(2:(ntaxa+1), function(k){ 
  l = strsplit(mashdist[[k]], split='\t')[[1]]
  if (k==2){ 
    t = c() 
  }else{ 
    t = as.numeric(l[2:length(l)]) 
  }
  return( c( t, rep(NA, (ntaxa+2-k)))) 
}))

assaccm = sapply(2:(ntaxa+1), function(k){
  nfseq = strsplit(mashdist[[k]], split='\t')[[1]][1]
  bnseq = basename(nfseq)
  nseq = gsub("(GC[AF]_[0-9]*\\.[1-9])_.*_genomic.fna.gz", "\\1", bnseq)
})

colnames(ml) = rownames(ml) = assaccm


print(nfassmetadata)
assmetadata = read.table(nfassmetadata, h=T, sep='\t', comment.char="")

longnames = sapply(assaccm, combineNames, metadata=assmetadata)

colnames(ml) = rownames(ml) = longnames
write.table(ml, sprintf('%s.tab', nfmashtriangle), col.names=T, row.names=T)

md = as.dist(ml)
sml = as.matrix(md)
mdd = hclust(md)
pmd = as.phylo(mdd)
pdf(sprintf('%s.hclust.pdf', nfmashtriangle), wid=30, hei=30)
plot(mdd)
dev.off()

write.tree(pmd, sprintf('%s.nwk', nfmashtriangle))

plottreeheatmap(pmd, sml, sprintf('%s.heatmap.pdf', nfmashtriangle))
if (!is.null(reftree)){
  t = reftree
  if (!is.null(treelab2ass)){
    assacct = sapply(reftree$tip.label, function(lab){
      i = which(treelab2ass[,2]==lab)
	  stopifnot(length(i)==1) #, sprintf("non-unique match to '%s' in metadata file", ass)
      return(treelab2ass[i,1])  
    })
    reftree$tip.label = longnames[assacct]
  }
  plottreeheatmap(reftree, sml, sprintf('%s_vs_%s.heatmap.pdf', nftree, basename(nfmashtriangle)))
}
