#################################################################################
### function to test module robustness
#	data - original data matrix/df
#		   if clustering animals based on behavior data then columns=samples, rows=behaviors
#
#	modules - numeric or character vector
#			  entries correspond to module assignments for each column of data
# robustness for one network!!
.bootstrapChiSqModuleOverlap = function (data, modules, nRuns=100, plot=T, ...) {
	chisqBoot = vector(mode='numeric', length=nRuns);
	nObserved = nrow(data);
	nModules = NULL # TODO rem
	tinyStatModules = NULL
	tinyStatRows = NULL
	normalStatModules = NULL
	normalStatRows = NULL
	for (run in 1:nRuns) {
		indexBoot = sample(1:nObserved, 3, replace=F);
		indexBoot = c(indexBoot, sample(1:nObserved, nObserved - 3, replace=T));
		pData = data[indexBoot, ];
		while(0 %in% apply(pData, 2, sum)) {
			indexBoot = sample(1:nObserved, 3, replace=F);
			indexBoot = c(indexBoot, sample(1:nObserved, nObserved - 3, replace=T));
			pData = data[indexBoot, ];
		}
		
		# run your usual module definition procedure on pData to get pModules
		# e.g.
	#	print(pData);
		pModules = cutreeDynamic(hclust(as.dist(1-cor(pData, use = 'p'))), distM = 1-cor(pData, use = 'p'), minClusterSize = 3, cutHeight = .995, verbose = 0);
		 if(length(table(pModules)) > 1) pModules = matchLabels(pModules, modules, pThreshold = .5); #NEW
		# tmpModules = if(length(table(pModules)) > 1) matchLabels(modules, pModules) else modules; #NEW
		nModules = c(nModules, length(table(pModules))); # TODO rem
		# you might find it useful to return all iterations of pModules in a matrix
		
		chisqBoot[run] = chisq.test(table(pModules, modules))$statistic;
		if(chisq.test(table(pModules, modules))$statistic < .1) {
			tinyStatModules = rbind(tinyStatModules, t(pModules));
			tinyStatRows = rbind(tinyStatRows, dimnames(pData)[[1]])
		} else {
			normalStatModules = rbind(normalStatModules, t(pModules));
			normalStatRows = rbind(normalStatRows, dimnames(pData)[[1]])
		}
		
		# depending on how long each run takes you may want to write chisqBoot to a file every n runs
	}
	
	if (plot) {
		titleMain = paste('mean=', signif(mean(chisqBoot), 5),
						  ', sd=', signif(sd(chisqBoot), 4),
						  ', cv=', signif(sd(chisqBoot)/mean(chisqBoot), 2),
						  sep=''
						  );
		sdLines = c((mean(chisqBoot)+sd(chisqBoot)), (mean(chisqBoot)-sd(chisqBoot)));
		hist(chisqBoot, main=titleMain, ...);
		abline(v=mean(chisqBoot), col='red');
		abline(v=sdLines, col='red', lty='dashed');
	}
	return(list(nModules = nModules, tinyStatModules = tinyStatModules, tinyStatRows = tinyStatRows, normalStatModules = normalStatModules, normalStatRows = normalStatRows)) # TODO rem
	return(chisqBoot);
}

#################################################################################

### function to plot heatmap of module overlaps across two networks
### depends on WGCNA::overlapTable and WGCNA::labeledHeatmap
#   modules1/modules2 - numeric or character vectors
#	thresh - logical indicating whether to only show numbers in cells that pass Bonferroni 
#   returnTable - logical indicating whether to return 
#   other args are hard-coded plot parameters that should work ok for networks with 20-50 modules
# compare 2 networks!!!
.plotModuleOverlaps = function(modules1, modules2, fcex=1.00, pcex=.7, fcex1=.7, pcex1=1.00, thresh=T, gamma=1, returnTable=T, ...) {
	overlap = overlapTable(modules1, modules2);
	
	numMat = -log10(overlap$pTable);
	numMat[numMat > 50] = 50;
	 
	textMat = paste(overlap$countTable, '\n', signif(overlap$pTable, 2));
	dim(textMat) = dim(numMat);
	 
	xlabels = paste('M', sort(unique(modules2)));
	ylabels = paste('M', sort(unique(modules1)));
	xSymbols = paste(sort(unique(modules2)), ': ', table(modules2), sep = '');
	ySymbols = paste(sort(unique(modules1)), ': ', table(modules1), sep = '');
	
	if (thresh) {
		thresh=.05/(ncol(textMat)*nrow(textMat));
		textMat[overlap$pTable>thresh] = ''
	}

	sizeGrWindow(7, 7); fp = FALSE
	par(mar = c(6, 7, 2, 1.0));
	labeledHeatmap(Matrix = numMat,
				   xLabels = xlabels, xSymbols = xSymbols,
				   yLabels = ylabels, ySymbols = ySymbols,
				   colorLabels = T, 
				   colors = blueWhiteRed(200,gamma=gamma)[80:200],
				   textMatrix = textMat, cex.text = pcex, setStdMargins = F,
				   cex.lab = fcex1,
				   xColorWidth = 0.08,
				   main = '', ...
				   );
				   
	if (returnTable) { return(overlap) }
}