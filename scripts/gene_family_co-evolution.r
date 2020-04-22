#/usr/bin/env Rscript --vanilla
library(parallel)

maxfreq = 1000
ncpus = 8

cargs = commandArgs(trailingOnly=T)

lnffamevents = readLines(cargs[1])
fams = sapply(lnffamevents, function(p){ strplit(basename(p), split='_samples.nhx')[[1]][1] })
F = length(fams)

famsizes = read.table(cargs[2], sep='\t', h=F, col.names=c('fam', 'size'))

lfameventtabs = lapply(lnffamevents, read.table, sep='\t', h=F, col.names=c('type', 'location', 'donor', 'freq'))
names(lfameventtabs) = fams
stopifnot(all(sort(fams)==sort(famisizes$fam)))

# truncate aggregated counts of multiple gene tree events with the same species tree address
for (fameventtab in lfameventtabs){
	fameventtab$freq[fameventtab$freq > maxfreq] = maxfreq
}

rowids = 1:(F-1)
lfamcoevol = mclapply(rowids, function(a){
	colids = (a+1):min(a+1, F)
	partialrow = sapply(colids, function(b){
		cotab = merge(lfameventtabs[[a]], lfameventtabs[[b]], by=1:3)
		cotab$coev = sqrt(cotab$freq.x * cotab$freq.y) / maxfreq
		# sum and scale by gene tree size
		Gtreesizes = (famsizes$size[famisizes$fam %in% fams[c(a, b)]] * 2) - 2
		return( sum(cotab$coev) / min(Gtreesizes) )
	})
	names(partialrow) = fams[colids]
	return(partialrow)
}, mc.cores=ncpus, mc.preschedule=F)
names(lfamcoevol) = fams[rowids]

mfamcoevol = matrix(NA, F, F)
for (a in 1:F){
 for (b in 1:F){
  if (a==b){ m = NA 
  }else{
   m = lfamcoevol[[fam[a]]][[fam[b]]]
   if (is.null(m)){ m = lfamcoevol[[fam[b]]][[fam[a]]] }
   mfamcoevol[a,b] = mfamcoevol[b,a] = m
  }
 }
}

names(mfamcoevol) = list(fams, fams)

heatmap(mfamcoevol)