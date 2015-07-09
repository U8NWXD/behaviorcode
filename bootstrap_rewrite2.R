bootstrap2independent_tmp = function(group1, group2, groupNames = c('group1', 'group2'), dataDescriptor = NULL,
									 trials = 10000, Func = 'mean', replace = T,
									 printResults = F, verbose = F, plots = T, ...) {
	data = .validateAndCleanBootstrapData(group1, group2, groupNames);	
	stat = eval(call(Func, data[[1]])) - eval(call(Func, data[[2]]));
	statsNULL = .buildNullDistribution(data, trials, Func, replace, dataDescriptor, verbose);
	midNULL = mean(statsNULL); # TODO add this and reflect to output
	reflect = midNULL - (stat - midNULL);
	p = .computePValue(stat, statsNULL, reflect, trials);
	
	if (plots) .makeBootstrapPlot(data, dataDescriptor, Func, stat, reflect, statsNULL, replace, trials, p, ...);

	output = list(stat = stat, p = p, null.dist = statsNULL, replacement = replace, data = data, parameters = match.call());
	if (printResults) {
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist and $parameters\n');
	}
	return(output);
}

.validateAndCleanBootstrapData = function(group1, group2, groupNames) {
	if (!is.numeric(c(group1,group2))) stop('DATA CONTAINS NON-NUMERIC VALUES.');
	if (sum(is.na(c(group1,group2))) > 0)
	{
		group1 = group1[!is.na(group1)];
		group2 = group2[!is.na(group2)];
		warning('NAs removed from one or both groups, check your data');
	}
	# TODO check if both groups are all-zero
	data = list(group1=group1, group2=group2);
	names(data) <- groupNames; # TODO check groupNames ok
	return(data);
}

.buildNullDistribution = function(data, trials, Func, replace, dataDescriptor, verbose) {
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
	statsNULL = numeric(length = trials); # try making it NA to test they all get replaced TODO
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
		warning("Statistic == 0; this may indicate a problem...", immediate. = T)
		p = 1;
	}
	return(p);
}

.makeBootstrapPlot = function(data, dataDescriptor, Func, stat, reflect, statsNULL, replace, trials, p,
							  col = 'grey', border = 'darkgrey', col.line = 'red', 
							  jitter = .15, pch = 21, boxLineMedian = F, outfile = NULL,
							  widthInInches = 10, heightInInches = 10) {
	if (!is.null(outfile)) jpeg(filename = outfile, width = widthInInches, heightInInches = 10, units = "in", quality = 100, res = 150, type = "quartz");
	if (is.null(dataDescriptor) || !is.character(dataDescriptor)) {
		dataDescriptor = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
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




