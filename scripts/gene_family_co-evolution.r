#!/usr/bin/env Rscript
library('parallel')
library('RColorBrewer')

midrange = function(x, ...){ mean(c(min(x, na.rm=T), max(x, na.rm=T)), ...) }

sminmedmax = function(mat){
  sapply(c('min', 'midrange', 'max'), function(fun){
	  sprintf("%.2g", get(fun)(mat, na.rm=T))
  })
}

maxfreq = 1000
evprobthresh = 0.5
ncpus = detectCores()

cargs = commandArgs(trailingOnly=T)

nflnffamevents = cargs[1]
nffamsizes = cargs[2]
nftestfams = cargs[3]
outprefix = cargs[4]
if (length(cargs)>4){
	only.eventtypes = strsplit(cargs[5], split=',')[[1]] 
}else{
	only.eventtypes = c('D', 'T', 'S', 'L')
}
outtag = paste(only.eventtypes, collapse='')
tmpdir = paste(outprefix, outtag, 'tmp_files')
hitbyquerydir = file.path(tmpdir, 'coevolscores_by_query_fam')
dir.create(hitbyquerydir, showWarnings=F, recursive=T)

lnffamevents = readLines(nflnffamevents)
fams = sapply(lnffamevents, function(p){ strsplit(basename(p), split='_samples')[[1]][1] })
f = length(fams)

famsizes = read.table(nffamsizes, sep='\t', h=F, col.names=c('fam', 'size'))

testfams = readLines(nftestfams)

lfameventtabs = lapply(lnffamevents, function(nffamevents){
	famevents = read.table(nffamevents, sep='\t', h=F, col.names=c('type', 'location', 'donor', 'freq'))
	famevents[famevents$type %in% only.eventtypes,]
})
names(lfameventtabs) = fams

# check data consistency
a = sapply(fams, function(fam){ fam %in% famsizes$fam })
stopifnot(all(a))
b = sapply(testfams, function(fam){ fam %in% famsizes$fam })
stopifnot(all(b))
c = sapply(testfams, function(fam){ fam %in% fams })
stopifnot(all(c))

controlfams = setdiff(fams, testfams)

# truncate aggregated counts of multiple gene tree events with the same species tree address
for (fameventtab in lfameventtabs){
	fameventtab$freq[fameventtab$freq > maxfreq] = maxfreq
}

rowids = 1:(f-1)
lfamcoevol = mclapply(rowids, function(a){
	colids = (a+1):f
	partialrow = sapply(colids, function(b){
#		cat(sprintf("\r# %d-%d (%s vs. %s)\t\n", a, b, fams[a], fams[b]))
		cotab = merge(lfameventtabs[[a]], lfameventtabs[[b]], by=1:3)
		cotab$coev = sqrt(cotab$freq.x * cotab$freq.y) / maxfreq
		# sum and scale by gene tree size
		Gtreesizes = (famsizes$size[famsizes$fam %in% fams[c(a, b)]] * 2) - 2
		return( sum(cotab$coev) / max(Gtreesizes) )
	})
	names(partialrow) = fams[colids]
	# save intermediary results
	write.table(partialrow, file=file.path(hitbyquerydir, paste(fams[colids], 'coevol-scores.tab'), sep='_'), quote=F, col.names=F, sep='\t')
	return(partialrow)
}, mc.cores=ncpus, mc.preschedule=F)
names(lfamcoevol) = fams[rowids]

lhighprobcoevt = mclapply(rowids, function(a){
	colids = (a+1):f
	partialrow = lapply(colids, function(b){
#		cat(sprintf("\r# %d-%d (%s vs. %s)\t\n", a, b, fams[a], fams[b]))
		cotab = merge(lfameventtabs[[a]], lfameventtabs[[b]], by=1:3)
		cotab$coev = sqrt(cotab$freq.x * cotab$freq.y) / maxfreq
		# sum and scale by gene tree size2
		return( cotab[cotab$coev>evprobthresh,] )
	})
	names(partialrow) = fams[colids]
	return(partialrow)
}, mc.cores=ncpus, mc.preschedule=F)
names(lhighprobcoevt) = fams[rowids]
write(capture.output(print(lhighprobcoevt[testfams])), file=paste(outprefix, outtag, 'high_probability_coevents', sep='.'))

mfamcoevol = matrix(NA, f, f)
for (a in 1:f){
 for (b in 1:f){
  if (a==b){ m = NA 
  }else{ if (a<b){ 
   m = lfamcoevol[[fams[a]]][[fams[b]]]
  }else{
   m = lfamcoevol[[fams[b]]][[fams[a]]]
  }}
  mfamcoevol[a,b] = mfamcoevol[b,a] = m
 }
}
colnames(mfamcoevol) = rownames(mfamcoevol) = fams
write.table(mfamcoevol, file=paste(outprefix, outtag, 'family_coevol_scores.mat', sep='.'), sep='\t')
mtestfamcoevol = mfamcoevol[testfams, testfams]
minscore = min(mfamcoevol)
for (i in 1:ncol(mfamcoevol)){ mfamcoevol[i,i] = minscore } # to keep colour scaling

pdf(paste(outprefix, outtag, 'family_coevol_scores.pdf', sep='.'), height=30, width=30)
heatmap(mfamcoevol, main='coevolution scores',
		col=colorRampPalette(brewer.pal(8, "Greens"))(25), scale='none')
legend(x="topleft", legend=sminmedmax(mfamcoevol), 
     fill=colorRampPalette(brewer.pal(8, "Greens"))(3))
dev.off()

write.table(mfamcoevol[testfams, testfams], file=paste(outprefix, outtag, 'test_family_coevol_scores.mat', sep='.'), sep='\t')

matscorepval = sapply(testfams, function(fam1){
  sapply(testfams, function(fam2){
	if (fam1==fam2){ s = 1 
	}else{ s = mfamcoevol[fam1, fam2] }
	bg1 = mfamcoevol[fam1, controlfams]
	bg2 = mfamcoevol[fam2, controlfams]
	bg = sort(c(bg1, bg2), decreasing=T)
	bg = bg[!is.na(bg)]
	for (k in 1:length(bg)){ if (s > bg[k]){ break } }
	return(k/length(bg))
  })
})
rownames(matscorepval) = colnames(matscorepval) = testfams
write.table(matscorepval, file=paste(outprefix, outtag, 'test_family_coevol_pvalues.mat', sep='.'), sep='\t')


pdf(paste(outprefix, outtag, 'test_family_coevol_scores_pvalues.pdf', sep='.'), height=15, width=15)
h = heatmap(mtestfamcoevol,
	 main='coevolution scores (test gene families only)',
	 col=colorRampPalette(brewer.pal(8, "Greens"))(25), scale='none', keep.dendro=T)
legend(x="topleft", legend=sminmedmax(mtestfamcoevol), 
     fill=colorRampPalette(brewer.pal(8, "Greens"))(3))
heatmap(matscorepval, Rowv=h$Rowv, Colv=h$Colv,
	 main='coevolution rank p-values',
	 col=colorRampPalette(brewer.pal(8, "Reds"))(25), scale='none')
legend(x="topleft", legend=sminmedmax(matscorepval), 
     fill=colorRampPalette(brewer.pal(8, "Reds"))(3))
dev.off()
	