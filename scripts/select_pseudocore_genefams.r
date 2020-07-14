#!/usr/bin/Rscript

minfracNgenomeshow=0.75

readInput = function(prompt){
	if (interactive() ){
		val <- readline(prompt=prompt)
	}else{
		cat(prompt);
		val <- readLines("stdin",n=1)
	}
	return(val)
}

promptselectmingenomes = function(P=-1, silent=F){
	cat("Please confirm or enter a new value for P [or stop by choosing 0]: \n")
	while (P < 0){
		P = as.numeric(readInput(prompt=ifelse(silent, "", "Please enter non-negative integer value for P: ")))
	}
	return(P)
}

getpseudocorefams = function(P, countsbyfam){
	pseudocorefams = names(countsbyfam[countsbyfam >= P])
	cat(sprintf("results in a set of %d pseudo-core unicopy gene families\n", length(pseudocorefams)))
	return(pseudocorefams)
}

selectMinGenomes = function(countmatrix, dirout, pseudocoremingenomes=-1, ngenomes=NA, interactive.X=FALSE, plot.PCoA=FALSE){
	countsbyfam = apply(countmatrix, 1, sum)
	if (is.numeric(ngenomes)){   N = ngenomes 
	}else{                       N = ncol(countmatrix) }
	cat(sprintf("number of unicopy gene families present in at least n genomes (out of %d):\n", N))
	cscbf = cumsum(rev(table(countsbyfam)))
	print(cscbf[intersect(names(cscbf), as.character(floor(N*minfracNgenomeshow):N))])
	if (is.null(pseudocoremingenomes)){ print("Select a value for P, the minimum number of genomes to be represented in pseudo-core unicopy gene families") }
	
	nloop = 0
	if (is.null(pseudocoremingenomes)){
		# interactive prompt
		p = promptselectmingenomes()
	}else{
		nloop = nloop + 1
		p = pseudocoremingenomes[nloop]
		if ( p <= 0 ) stop("needs a non-null value for 'pseudocoremingenomes'")
		print(sprintf("test value %d for P, the minimum number of genomes to be represented in pseudo-core unicopy gene families", p), quote=F)
	}
	pseudocorefams = NULL
	pseudocore = list()
	pcmg = -1
	if (interactive.X){ X11(width=16, height=10) }
	while (p!=pcmg | nloop<=length(pseudocoremingenomes)){
		if (p<N){ pseudocorerad = sprintf("pseudo-core-%d-unicopy", p)
		}else{  pseudocorerad = "strict-core-unicopy" }
		nftabout = file.path(dirout, sprintf("%s_families.tab", pseudocorerad))
		nfpdfout = file.path(dirout, sprintf("%s_families.pdf", pseudocorerad))
		if (!interactive.X){ pdf(nfpdfout, width=40, height=40) }
		pcmg = p
		pseudocorefams = getpseudocorefams(p, countsbyfam)
		if (length(pseudocorefams)>0){
			cat("plotting heatmap... ")
			pseudocoremat = countmatrix[pseudocorefams,]
			heatmap(t(pseudocoremat), breaks=c(-0.5, 0.5, 1.5, N), col=c('white', 'black', 'red'), scale='none')
			if (interactive.X){
				cat("please verify that the distribution of markers per species is not too skewed (counts per species, white: 0, black: 1, red: >1)\n")
				message("Press Return To Continue") ; invisible(readLines("stdin", n=1))
			}
			if (plot.PCoA){
				cat("computing and plotting PCoA of gene species based on presence/absence... ")
				count.coa = ade4::dudi.coa(pseudocoremat, scannf=F, nf=2)
				if (interactive.X){ message("Press Return To Continue") ; invisible(readLines("stdin", n=1)) }
				s.label(count.coa$c1)
			}
			nmissing = apply(pseudocoremat, 2, function(x){ length(which(!x)) })
			barplot(nmissing, horiz=T, las=2, ylab='Species label', xlab='Nb. missing gene markers')
			if (interactive.X){ message("Press Return To Continue") ; invisible(readLines("stdin", n=1)) }
			barplot(nmissing[order(nmissing, decreasing=T)[1:min(20, N)]], horiz=T, las=2, ylab='Species label', xlab='Nb. missing gene markers')
			if (!interactive.X){ dev.off() }
		}
		write(pseudocorefams, file=nftabout)
		cat(sprintf("Written list of %d pseudo-core unicopy gene families (with min. genome nb. = %d) and graphical representation of their distribution at:\n%s\n%s\n",
		 length(pseudocorefams), p, nfpdfout, nftabout))
		pseudocore[[as.character(p)]] = pseudocorefams
		nloop = nloop + 1
		if (nloop < length(pseudocoremingenomes)){
			p = pseudocoremingenomes[nloop]
		}else{ if (is.null(pseudocoremingenomes)){
			# interactive prompt
			p = promptselectmingenomes()
		}}
		# else just keep value of p, which will lead to p == pcmg and breaking the while loop
		if (p==0){ break }
	}
	if (interactive.X){ dev.off() }
	cat(sprintf("Selected %d as value of P\n", p))
	return(pseudocore)
}

cargs = commandArgs(trailingOnly=TRUE)

nffamgenomemat = cargs[1]
nflasscode = cargs[2]
dirout = cargs[3]
if (length(cargs) > 3){
	pseudocoremingenomes = as.numeric(cargs[4])
	if (!is.na(pseudocoremingenomes)){
		cat(sprintf("try value of pseudocoremingenomes = %d\n", pseudocoremingenomes))
	}else{
		if (cargs[4]==''){
			pseudocoremingenomes = NULL
			cat(sprintf("no input for pseudocoremingenomes; turns on INTERACTIVE mode\n", pseudocoremingenomes))
		}
		print(sprintf("extract values of pseudocoremingenomes from file '%s'", cargs[4]), quote=F)
		pseudocoremingenomes = as.numeric(readLines(cargs[4]))
		cat(sprintf("try values of pseudocoremingenomes = %s\n", paste(pseudocoremingenomes, collapse=', ')))
	}
}else{
	pseudocoremingenomes = NULL
}
if (length(cargs) > 4){
	interactive.X = as.logical(toupper(cargs[5]))
	if (is.na(interactive.X)){ interactive.X = as.logical(as.numeric(cargs[5])) }
}else{
	interactive.X = FALSE
}
if (length(cargs) > 5){
	nfrestrictlist = cargs[6]
}else{
	nfrestrictlist = NULL
}
if (length(cargs) > 6){
	restricttag = cargs[7]
}else{
	if (!is.null(nfrestrictlist)){
		restricttag = strsplit(basename(nfrestrictlist), split='_genome_codes')[[1]][1]
	}else{
		restricttag = NULL
	}
}
cat("Loading matrix of gene families counts in genomes...\n")
genocount = data.matrix(read.table(file=nffamgenomemat, sep='\t', comment.char='', header=T, row.names=1, check.names=T))
cat("Loading correspondence table of assembly accessions to genome codes...\n")
lasscode = read.table(nflasscode, row.names=1, stringsAsFactors=F, sep='\t', comment.char='')
# get names to be consistent
rownames(lasscode) = make.names(rownames(lasscode))
colnames(genocount) = lasscode[colnames(genocount),1]

if (!is.null(nfrestrictlist)){
	restrictgenomelist = readLines(nfrestrictlist)
	genocount = genocount[,restrictgenomelist]
    cat("Saving restricted matrix of gene families counts in genomes...\n")
	write.table(genocount, file=sprintf("%s_restricted_%s", nffamgenomemat, restricttag), sep='\t', quote=F)
}

onlyunicopy = apply(genocount, 1, function(x){ max(x)==1 })
genocountunicopy = genocount[onlyunicopy,]

pseudocore = selectMinGenomes(genocountunicopy, dirout, pseudocoremingenomes=pseudocoremingenomes, interactive.X=interactive.X)
nfdataout = file.path(dirout, "pseudo-core-all.RData")
save(genocountunicopy, pseudocore, file=nfdataout)
cat(sprintf("Saved data in file: '%s'.\n", nfdataout))

P = names(pseudocore)[length(pseudocore)]
cat(sprintf("Final choice of %d pseudo-core unicopy gene families (present in at least %d genomes).\n", length(pseudocore[[P]]), as.numeric(P)))

write(sprintf("pseudocoremingenomes=%d", as.numeric(P)), stderr())
