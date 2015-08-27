# The functions in this file are heavily based on the code written by Austin in
# "bootstrap_tests_June2013_STABLE.R". A few changes were made related to moving
# things out of for-loops to improve efficiency, and functions were decomposed for
# readability.



# Runs a bootstrap test on two numeric data vectors (<group1> and <group2>) representing
# independent groups (names should be given in <groupNames>). Uses <Func> (either "mean"
# or "median") to compute the stats.
#
# If <plots> is TRUE, a plot is output to the active Quartz window, or a jpeg if <outfile>
# is specified. Graphical parameters that can be passed to the plotting function include:
# <col> (default 'grey'), the color of the histogram bars;
# <border> (default 'darkgrey'), the color of plot borders;
# <col.line> (default 'red'), the color of the lines showing means or medians, and <stat>
#   on the plot of the bootstrap results;
# <jitter> (default .15), the amount of jitter that the scatterplot of data on top of the boxplot has;
# <pch> (default 21, an open circle), the type of point used in this scatterplot;
# <boxLineMedian> (default FALSE), should a line be drawn in the box plot at the median? FALSE makes
#   it at the mean instead;
# <widthInInches> (default 10) and <heightInInches> (default 10), the dimensions of the output jpeg.
bootstrap2independent = function(group1, group2, groupNames = c('group1', 'group2'), dataDescriptor = NULL,
								 trials = 10000, Func = 'mean', replace = T,
								 printResults = F, verbose = F, plots = T, ...) {
	data = .validateAndCleanBootstrapData(group1, group2, groupNames, paired = F);
	if (is.null(data)) return(list(p = NA));	
	stat = eval(call(Func, data[[1]])) - eval(call(Func, data[[2]]));
	statsNULL = .buildNullDistributionIndependent(data, trials, Func, replace, dataDescriptor, verbose);
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	p = .computePValue(stat, statsNULL, reflect, trials);
	
	if (plots) .makeBootstrapPlotIndependent(data, dataDescriptor, Func, stat, reflect, statsNULL, replace, trials, p, ...);

	output = list(stat = stat, stat.reflect = reflect, p = p, null.dist = statsNULL, null.dist.mean = midNULL,
				  replacement = replace, data = data, parameters = match.call());
	if (printResults) {
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist and $parameters\n');
	}
	return(output);
}

# Runs a bootstrap test on two numeric data vectors (<condition1> and <condition2>) representing
# measurements taken at two timepoints on the same group of subjects. (names should be given in
# <conditionNames>). Uses <Func> (either "mean" or "median") to compute the stats.
#
# If <plots> is TRUE, a plot is output to the active Quartz window, or a jpeg if <outfile>
# is specified. Graphical parameters that can be passed to the plotting function include:
# <col> (default 'grey'), the color of the histogram bars;
# <border> (default 'darkgrey'), the color of plot borders;
# <col.line> (default 'red'), the color of the lines showing means or medians, and <stat>
#   on the plot of the bootstrap results;
# <pch> (default 21, an open circle), the type of point used in this scatterplot;
# <cex.lab> and <cex.axis> (default 1.2), the magnifications of labels and axes;
# <widthInInches> (default 12) and <heightInInches> (default 5), the dimensions of the output jpeg.
bootstrap2paired = function(condition1, condition2, conditionNames = c('condition1', 'condition2'), dataDescriptor = NULL,
							trials = 10000, Func = 'mean', abs.diffs = F,
							printResults = F, verbose = F, plots = T, ...) {
	data = .validateAndCleanBootstrapData(condition1, condition2, conditionNames, paired = T);
	if (is.null(data)) return(list(p = NA));
	
	diffs = data[[1]] - data[[2]];
	stat = eval(call(Func, diffs));
	if (abs.diffs) stat = abs(stat);
	
	statsNULL = .buildNullDistributionPaired(diffs, trials, Func, abs.diffs, dataDescriptor, verbose);
	if (abs.diffs) {
		p = sum(statsNULL > stat) / trials;
	} else {
		midNULL = mean(statsNULL);
		reflect = midNULL - (stat - midNULL);
		p = .computePValue(stat, statsNULL, reflect, trials);
	}
	
	if (plots) .makeBootstrapPlotPaired(data, diffs, dataDescriptor, Func, stat, reflect, statsNULL, abs.diffs, trials, p, ...)
	
	output = list(stat = stat, p=p, null.dist = statsNULL, data = data.frame(data, diffs = diffs), parameters = match.call());
	if (!abs.diffs) {
		output$stat.reflect <- reflect;
		output$null.dist.mean <- midNULL;
	}
	if (printResults) {
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist and $parameters\n');
	}
	return(output);
}

# Helper function for bootstrap2independent() and bootstrap2paired()
# Checks the data given in <group1> and <group2> to ensure it is numeric, nonempty, and not all 0s or NAs.
# If <paired>, the data is also checked to be sure it is the same length.
# Either way, NAs are removed, groupNames is checked to be a character vector of length 2, and the cleaned
# data is returned.
# If problems are found, this function may throw an error, throw a warning, or cause the main function
# to return early with p = NA.
.validateAndCleanBootstrapData = function(group1, group2, groupNames, paired) {
	if (!is.numeric(c(group1,group2))) stop('DATA CONTAINS NON-NUMERIC VALUES.');
	if (paired && length(group1) != length(group2)) stop('GROUPS ARE DIFFERENT SIZES...\n   THIS DATA SHOULD BE PAIRED!!!');
	if (sum(is.na(c(group1,group2))) > 0) {
		NAcheck1 = is.na(group1);
		NAcheck2 = is.na(group2);
		if (paired) NAcheck1 = NAcheck2 = NAcheck1 | NAcheck2
		group1 = group1[!NAcheck1];
		group2 = group2[!NAcheck2];
		if (paired) warning('NAs in one/both conditions, corresponding data points removed from BOTH')
		else warning('NAs removed from one or both groups, check your data');
	}
	
	if (length(group1) == 0 || length(group2) == 0) {
		warning('All data was NA in one group or the other. Aborting.', immediate. = T);
		return(NULL);
	} else if (sum(group1) + sum(group2) == 0) {
		warning('All data in both groups is 0. Aborting.', immediate. = T);
		return(NULL);
	}
	data = list(group1=group1, group2=group2);
	
	if (!is.character(groupNames)) {
		warning("groupNames is not a character vector. Coercing to character...", immediate.=T);
		groupNames = as.character(groupNames);
	}
	if (length(groupNames) < 2) {
		stop("groupNames must be length 2.")
	} else if (length(groupNames) > 2) {
		warning('groupNames has more than two elements. Using the first two names ("', groupNames[1], '" and "', groupNames[2], '")...', immediate.=T);
		groupNames = groupNames[1:2];
	}
	names(data) <- groupNames;
	
	return(data);
}

# Helper function for bootstrap2independent()
# Builds the null distribution for independent data.
.buildNullDistributionIndependent = function(data, trials, Func, replace, dataDescriptor, verbose) {
	if (verbose) {
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: difference of group ', Func, 's\n', sep = '');
		cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
	}
	
	boxNULL = c(data[[1]], data[[2]]);
	group1Length = length(data[[1]]);
	group2Length = length(data[[2]]);
	fxn = get(Func);
	statsNULL = numeric(length = trials); # this improves speed rather than expanding a vector 10000 times
	for (trial in 1:trials) {
		if (verbose && trial %% 1000 == 0) {
			if (trials <= 20000) cat('  Run ', trial, '\n', sep = '')
			else if (trial %% 10000 == 0) cat('  Run ', trial, '\n', sep = '');
		}
		pseudo_group1 = sample(boxNULL, group1Length, replace = replace);
		pseudo_group2 = sample(boxNULL, group2Length, replace = replace);
		pseudo_stat_null = fxn(pseudo_group1) - fxn(pseudo_group2);
		statsNULL[trial] = pseudo_stat_null;
	}
	
	return(statsNULL);
}

# Helper function for bootstrap2paired()
# Builds the null distribution for paired data.
.buildNullDistributionPaired = function(diffs, trials, Func, abs.diffs, dataDescriptor, verbose) {
	if (verbose) {
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: ', Func, ' of individual differences\n\n', sep = '');
	}
	
	fxn = get(Func);
	numDiffs = length(diffs);
	statsNULL = numeric(length = trials); # this improves speed rather than expanding a vector 10000 times
	for (trial in 1:trials) {
		if (verbose && trial %% 1000 == 0) {
			if (trials <= 20000) cat('  Run ', trial, '\n', sep = '')
			else if (trial %% 10000 == 0) cat('  Run ', trial, '\n', sep = '');
		}
		signs = sample(c(-1,1), numDiffs, replace = T)
		pseudo_diffs = signs * diffs;
		pseudo_stat_null = fxn(pseudo_diffs);
		statsNULL[trial] = pseudo_stat_null;
	}
	if (abs.diffs) statsNULL = abs(statsNULL);
	return(statsNULL);
}

# Helper function for bootstrap2independent() and bootstrap2paired()
# Computes the p-value from statsNULL.
.computePValue = function(stat, statsNULL, reflect, trials) {
	if (stat < 0) { 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} else if (stat > 0) {
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	} else {
		warning("Statistic == 0; this may indicate a problem...")
		p = 1;
	}
	return(p);
}

# Helper function for bootstrap2Independent().
# Makes a plot showing histograms of both groups, a histogram of the results of the bootstrap run,
# and box-and-whisker plots overlaid with individual data points showing the data from each group.
.makeBootstrapPlotIndependent = function(data, dataDescriptor, Func, stat, reflect, statsNULL, replace, trials, p,
										 col = 'grey', border = 'darkgrey', col.line = 'red', 
										 jitter = .15, pch = 21, boxLineMedian = F, outfile = NULL,
										 widthInInches = 10, heightInInches = 10) {
	if (!is.null(outfile)) jpeg(filename = outfile, width = widthInInches, height = heightInInches, units = "in", quality = 100, res = 150, type = "quartz");
	if (is.null(dataDescriptor) || !is.character(dataDescriptor)) {
		dataDescriptor = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone', immediate. = T);
	}
	names(data) <- paste(names(data), ' (n=', unlist(lapply(data, length)),')', sep = '');
	all_data = c(data[[1]], data[[2]]);
	
	par(mfrow = c(2,2));

	# make histograms
	hist_group1 = hist(data[[1]], breaks = length(all_data) / 2, plot = F);
	hist_group2 = hist(data[[2]], breaks = length(all_data) / 2, plot = F);
	ymax = max(c(hist_group1$counts, hist_group2$counts));
	xmin = min(c(hist_group1$breaks, hist_group2$breaks));
	xmax = max(c(hist_group1$breaks, hist_group2$breaks));
	
	plot(hist_group1, main = names(data)[1], ylim = c(0, ymax), xlim = c(xmin, xmax), 
		 col = col, border = border, xlab = dataDescriptor);
	abline(v = eval(call(Func, data[[1]])), col = col.line);
	plot(hist_group2, main = names(data)[2], ylim = c(0, ymax), xlim = c(xmin, xmax), 
		 col = col, border = border, xlab = dataDescriptor);
	abline(v = eval(call(Func, data[[2]])), col = col.line);
	
	# plot histogram of null distribution
	hist(statsNULL, main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
		 xlab = paste('statistic (', Func, 's, replace = ', substr(replace, 1, 1), ')', sep = ''),
		 col = col, border = border, breaks = floor(trials / 100));
	abline(v = stat, col = col.line);
	abline(v = reflect, col = col.line, lty = 'dashed');
	
	# boxplots
	toPlot = all_data;
	grp = c(rep(names(data)[1], length(data[[1]])), rep(names(data)[2], length(data[[2]])));
	medlty = if(boxLineMedian) "solid" else "blank";
	toPaste = if(p == 0) 'p < 1e-5' else paste('p = ', signif(p, 2), sep = '');
	boxplot(toPlot ~ grp, ylab = dataDescriptor, medlty = medlty, main = toPaste);
	if (!boxLineMedian) {
		segments(0.6, mean(data[[1]]), 1.4, mean(data[[1]]), lwd = 2);
		segments(1.6, mean(data[[2]]), 2.4, mean(data[[2]]), lwd = 2);
	}
	stripchart(toPlot ~ grp, vertical = T, add = T, method = 'jitter', jitter = jitter,
			   pch = pch, bg = col, cex = 1.5);
			   
	if (!is.null(outfile)) dev.off();
}

# Helper function for bootstrap2Paired().
# Makes a plot showing histograms of the individual differences between condition1 and condition2,
# a histogram of the results of the bootstrap run, and box-and-whisker plots overlaid with paired
# data points showing the data from each group.
.makeBootstrapPlotPaired = function(data, diffs, dataDescriptor, Func, stat, reflect, statsNULL, abs.diffs, trials, p,
									col = 'grey', border = 'darkgrey', col.line = 'red', 
									pch = 21, outfile = NULL, cex.lab = 1.2, cex.axis = 1.2,
									widthInInches = 12, heightInInches = 5, ...) {
	if (!is.null(outfile)) jpeg(filename = outfile, width = widthInInches, height = heightInInches, units = "in", quality = 100, res = 150, type = "quartz");
	if (is.null(dataDescriptor) || !is.character(dataDescriptor)) {
		dataDescriptor = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone', immediate. = T);
	}
	
	par(mfrow = c(1,3));

	# make histogram of diffs
	hist_diffs = hist(diffs, breaks = length(diffs), plot = F);
	xmax = max(abs(hist_diffs$breaks));
	plot(hist_diffs, main = paste('Individual differences with ', Func, '\n(n=', length(diffs), ')', sep = ''), xlim = c(-xmax, xmax), 
		 xlab = paste(dataDescriptor, ' (', names(data)[1], ' - ', names(data)[2], ')', sep = ''),
		 col = col, border = border, cex.lab = cex.lab, cex.axis = cex.axis);
	abline(v = eval(call(Func, diffs)), col = col.line);
	
	# plot histogram of null distribution
	xlabNULL = if (abs.diffs) paste('statistic: abs(', Func, ' of individual differences)', sep = '')
			   else paste('statistic: ', Func, ' of individual differences', sep = '');
	hist(statsNULL, main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
		 xlab = xlabNULL, col = col, border = border, breaks = floor(trials / 100),
		 cex.lab = cex.lab, cex.axis = cex.axis);
	abline(v = stat, col = col.line);
	if (!abs.diffs) abline(v = reflect, col = col.line, lty = 'dashed');
	
	# boxplots
	toPlot = as.data.frame(data);
	toPaste = if(p == 0) paste('p < 1e-5 (n=', length(diffs), ')', sep = '') else paste('p = ', signif(p, 2), ' (n=', length(diffs), ')', sep = '');
	boxplot(toPlot, ylab = dataDescriptor, main = toPaste, frame.plot = F, cex.lab = cex.lab, cex.axis = cex.axis, ...);
	stripchart(toPlot, vertical = T, add = T, pch = pch, bg = col, cex = 1.5);
	segments(1, toPlot[,1], 2, toPlot[,2], col = border);
			   
	if (!is.null(outfile)) dev.off();
}





