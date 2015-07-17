library(dynamicTreeCut)



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
	
		for (clustStyle in c("complete", "average")) {
			hc = hclust(as.dist(1 - correlationsMat), method = clustStyle);
			jpeg(paste(paste(outPref, name, "hclust", clustStyle, sep = "_"), "jpg", sep = '.'), width = 1024, height = 1024)
			plot(hc, main = paste(name, ' (', clustStyle, ')', sep = ''));
			dev.off();
			
			# for (minClusterSize in 3:(dim(correlationsMat)[1] / 3)) {
			for (minClusterSize in 3:5) {
				moduleNums = cutreeDynamic(hc, distM = 1 - correlationsMat, minClusterSize = minClusterSize, cutHeight = .995);
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