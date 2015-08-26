library(dynamicTreeCut)
library(WGCNA)

# logLabels should be c('ND', 'ND', 'D', 'D', 'ND', 'ASC', 'D', 'ND')
# or c('DOM1 (0 fish killed)', 'DOM2 (3 fish killed)', 'DOM2 (12 fish killed)') or the like. basically the labels for the dendro.
.clustersEtc = function(logList, clusterOn = 'counts', logLabels = NULL, cor.use = 'all.obs', hclust.method = 'complete', cutree.minClusterSize = NA, colors = NULL) {
	if (is.na(cutree.minClusterSize)) {
		cutree.minClusterSize = max(length(logList) / 6, 3);
	}
	if (is.null(colors)) colors = c('darkgrey', 'red', 'green', 'blue', 'orange', 'purple', 'darkgreen', 'gold')
	else colors = c('darkgrey', colors)
	
	
	if (clusterOn == 'counts') {
		dataMat = .extractBasicStats(logList, .behnames(logList), NULL)[[1]]
	} else if (clusterOn %in% c('tps', 'TPs', 'transprobs', 'transitional probabilities')) {
		dataMat = .makeTPMatrix(lapply(logList, function(d) {.getProbabilityMatrix(d$behavior, byTotal=F)}), .behnames(logList), FALSE)
	} else stop('THIS FUNCTION KNOWS NOTHING OF THESE "', clusterOn, '" YOU SPEAK OF');
	
	if (!is.null(logLabels)) colnames(dataMat) <- logLabels;
	
	corMat = cor(dataMat, use = cor.use);
	cluster = hclust(as.dist(1 - corMat), method = hclust.method);
	moduleNumbers = cutreeDynamic(cluster, minClusterSize = cutree.minClusterSize, distM = 1 - corMat, verbose = 0)
	if (max(moduleNumbers) < length(colors)) moduleColors = colors[moduleNumbers + 1]
	else {
		warning('Not enough colors. The names will be ugly now sorry.', immediate. = T)
		moduleColors = c('darkgrey', rainbow(max(moduleNumbers)))[moduleNumbers + 1];
	}
	plotDendroAndColors(cluster, moduleColors);
	return(list(data = dataMat, correlations = corMat, dendrogram = cluster, module.assignments = moduleNumbers, module.colors = moduleColors, call = match.call()))
}





.generateLotsOfDendrograms = function(dataMatrices, outPref) {
	output = list()
	for (i in 1:length(dataMatrices)) {
		dataMatrix = dataMatrices[[i]];
		name = names(dataMatrices)[i];
		moduleAssignments = NULL;
		correlationsMat = cor(dataMatrix, use = "p");
		jpeg(paste(outPref, name, "heatmap.jpg", sep = "_"), width = 1024, height = 1024)
		par(oma = c(20,0,3,15));
		heatmap(correlationsMat, symm = T, main = name);
		dev.off();
		
		print(name)
	
		for (clustStyle in c("complete", "average")) {
			hc = hclust(as.dist(1 - correlationsMat), method = clustStyle);
			jpeg(paste(paste(outPref, name, "hclust", clustStyle, sep = "_"), "jpg", sep = '.'), width = 1024, height = 1024)
			plot(hc, main = paste(name, ' (', clustStyle, ')', sep = ''));
			dev.off();
			
			# for (minClusterSize in 3:(dim(correlationsMat)[1] / 3)) {
			for (minClusterSize in 3) {
				moduleNums = cutreeDynamic(hc, distM = 1 - correlationsMat, minClusterSize = minClusterSize, cutHeight = .995, verbose = 0);
				moduleColors = c('darkgrey', rainbow(max(moduleNums)))[moduleNums + 1];
				jpeg(paste(paste(outPref, name, "modules", clustStyle, "minClusterSize", minClusterSize, sep = "_"), "jpg", sep = '.'), width = 1024, height = 1024)
				plotDendroAndColors(hc, moduleColors, main = paste(name, ' (', clustStyle, ', minClusterSize = ', minClusterSize, ')', sep = ''));
				dev.off();
				moduleAssignments = rbind(moduleAssignments, moduleNums);
				dimnames(moduleAssignments)[[1]][dim(moduleAssignments)[[1]]] <- paste(clustStyle, minClusterSize); 
			}
		}
		dimnames(moduleAssignments)[[2]] <- dimnames(correlationsMat)[[2]];
		output$.tmp = moduleAssignments;
		names(output)[which(names(output) == ".tmp")] <- name;
	}
	return(output);
}

.makeBoxplot = function(data, lognames, ...) {
	verboseBoxplot(data, gsub("/.*$", '', gsub("darkgrey", "zzz", lognames)), border = c('darkblue', 'gold', 'grey', 'transparent'), names = c('dark blue', 'gold', 'unassigned', 'not included'), xlab = 'Module', notch = F, frame.plot = F, cex.axis = 1.1, ...)
	pointsStaggered(1, data[gsub("/.*$", '', lognames) == "darkblue"], "darkblue")
	pointsStaggered(2, data[gsub("/.*$", '', lognames) == "gold"], "gold")
	pointsStaggered(3, data[gsub("/.*$", '', lognames) == "grey"], "grey")
	pointsStaggered(4, data[gsub("/.*$", '', lognames) == "darkgrey"], "grey35")
	# points(rep(1,length(data[gsub("/.*$", '', lognames) == "darkblue"])), data[gsub("/.*$", '', lognames) == "darkblue"], pch = 16, cex = 2, col = "darkblue")
	# points(rep(2,length(data[gsub("/.*$", '', lognames) == "gold"])), data[gsub("/.*$", '', lognames) == "gold"], pch = 16, cex = 2, col = "gold")
	# points(rep(3,length(data[gsub("/.*$", '', lognames) == "grey"])), data[gsub("/.*$", '', lognames) == "grey"], pch = 16, cex = 2, col = "grey")
	# points(rep(4,length(data[gsub("/.*$", '', lognames) == "darkgrey"])), data[gsub("/.*$", '', lognames) == "darkgrey"], pch = 16, cex = 2, col = "grey35")
}

pointsStaggered = function(x, y, color, pointsspace = .05) {
	freqtable = table(y);
	n = 1;
	npointsPlotted = 0
	while (sum(freqtable >= n)) {
		toPlot = as.numeric(names(freqtable)[freqtable == n])
		exes = x + (((1:n) - (n+1)/2) * pointsspace * 2);
		for (i in exes) {
			points(rep(i, length(toPlot)), toPlot, pch = 16, cex = 2, col = color);
			npointsPlotted = npointsPlotted + length(toPlot)
		}
		n = n+1;
	}
	print(npointsPlotted)
} 


par(mfrow = c(4,1))
par(oma = c(0, 1, 0, 0))
.makeBoxplot(x$durations[1,], names(counts), ylab = "Time attending to females (s)")
.makeBoxplot(x$durations[2,], names(counts), ylab = "Time attending to males (s)")
.makeBoxplot(x$durations[4,], names(counts), ylab = "Time maintaining territory (s)")
.makeBoxplot(counts, names(counts), ylab = "Total number of behaviors")
dev.off()