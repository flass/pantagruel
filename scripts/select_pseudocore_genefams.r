#!/usr/bin/Rscript --vanilla --silent
#~ library('ade4')

#~ options(width = 160)

minfracNgenomeshow=0.7

readInput = function(prompt){
	if (interactive() ){
		val <- readline(prompt=prompt)
	}else{
		cat(prompt);
		val <- readLines("stdin",n=1)
	}
	return(val)
}

selectmingenomes = function(pseudocoremingenomes=NA, silent=F){
	while (is.na(pseudocoremingenomes) || pseudocoremingenomes < 1){
		pseudocoremingenomes = as.numeric(readInput(prompt="Please enter non-null integer value for minimum of genomes represented in pseudo-core unicopy gene families: "))
	}
	if (!silent){ cat(sprintf("Selected %d as the minimum number of genomes to be represented in pseudo-core unicopy gene families\n", pseudocoremingenomes)) }
	return(pseudocoremingenomes)
}

getpseudocorefams = function(pseudocoremingenomes, countsbyfam){
	pseudocorefams = names(countsbyfam[countsbyfam >= pseudocoremingenomes])
	cat(sprintf("results in a set of %d pseudo-core unicopy gene families\n", length(pseudocorefams)))
	return(pseudocorefams)
}

selectMinGenomes = function(countmatrix, dirout, pseudocoremingenomes=NA, ngenomes=NULL){
	countsbyfam = apply(countmatrix, 1, sum)
	if (is.null(ngenomes)){ N = max(countsbyfam) }else{ N = ngenomes }
	cat(sprintf("number of unicopy gene families present in at least n genomes (out of %d):\n", N))
	print(cumsum(rev(table(countsbyfam)))[as.character(floor(N*minfracNgenomeshow):N)])

	pseudocoremingenomes = selectmingenomes(pseudocoremingenomes)	
	
	pcmg = -1
#~ 	X11(width=16, height=10)
	while (pseudocoremingenomes != pcmg){
		nftabout = file.path(dirout, sprintf("pseudo-core-%d-unicopy_families.tab", pseudocoremingenomes))
		nfpdfout = file.path(dirout, sprintf("pseudo-core-%d-unicopy_families.pdf", pseudocoremingenomes))
		pdf(nfpdfout, width=30, height=20)
		pcmg = pseudocoremingenomes
		pseudocorefams = getpseudocorefams(pseudocoremingenomes, countsbyfam)
		cat("please verify that the distribution of markers per species is not too skewed (counts per species, white: 0, black: 1, red: >1)\n")
		cat("plotting heatmap... ")
		heatmap(countmatrix[pseudocorefams,], breaks=c(-0.5, 0.5, 1.5, N), col=c('white', 'black', 'red'), scale='none')
		message("Press Return To Continue") ; invisible(readLines("stdin", n=1))
#~ 		cat("computing and plotting PCoA of gene species based on presence/absence... ")
#~ 		count.coa = dudi.coa(countmatrix[pseudocorefams,], scannf=F, nf=2)
#~ 		s.label(count.coa$c1)
#~ 		message("Press Return To Continue") ; invisible(readLines("stdin", n=1))
		nmissing = apply(countmatrix[pseudocorefams,], 2, function(x){ length(which(!x)) })
		barplot(nmissing, las=2, xlab='Species label', ylab='Nb. missing gene markers')
		message("Press Return To Continue") ; invisible(readLines("stdin", n=1))
		barplot(nmissing[order(nmissing, decreasing=T)[1:min(20, N)]], las=2, xlab='Species label', ylab='Nb. missing gene markers')
		cat("please confirm value for minimum of genomes represented in pseudo-core unicopy gene families: \n")
		pseudocoremingenomes = selectmingenomes(silent=T)
		dev.off()
		write(pseudocorefams, file=nftabout)
	}
	#~ 	exportcmd = sprintf("export pseudocoremingenomes=%d", pseudocoremingenomes)
	#~ system(exportcmd) ; cat(sprintf("system call: %s\n", exportcmd))
#~ 	dev.off()
	return(list(mingenomes=pseudocoremingenomes, fams=pseudocorefams))
}

#~ ngenomes = as.numeric(Sys.getenv('ngenomes'))
#~ pseudocoremingenomes = as.numeric(Sys.getenv('pseudocoremingenomes'))
#~ protali = Sys.getenv('protali')
#~ nflasscode =  file.path(Sys.getenv('database'), 'cds_codes.tab')
#~ dirout = protali
#~ nffamgenomemat = file.path(protali, 'full_families_genome_counts-noORFans.mat')

cargs = commandArgs()

ngenomes = cargs[1]
nffamgenomemat = cargs[2]
nflasscode = cargs[3]
dirout = cargs[4]
if (length(cargs) > 4){
	pseudocoremingenomes = cargs[5]
}else{
	pseudocoremingenomes = NA
}
cat("Loading matrix of gene families counts in genomes...\n")
genocount = data.matrix(read.table(file=nffamgenomemat))
onlyunicopy = apply(genocount, 1, function(x){ max(x)==1 })

pseudocore = selectMinGenomes(genocount[onlyunicopy,], dirout, pseudocoremingenomes=pseudocoremingenomes, ngenomes=ngenomes)


cat(sprintf("Written list of pseudo-core unicopy gene families (with min. genome nb. = %d) at '%s'.\n", pseudocore$mingenomes, nfout))

quit(status=pseudocoremingenomes, save='no')
