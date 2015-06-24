#############
####  ####
#############
confidenceInterval = function(data, 
							  interval = 0.95, 
							  trials = 10000, 
							  Func = 'mean', 
							  plots = T, 
							  col = 'grey', 
							  border = 'darkgrey', 
							  line.col = 'red',
							  pch = 21
							  )
{
	if (!is.numeric(data) | !is.vector(data))
	{
		stop('CHECK THAT DATA IS A NUMERIC VECTOR');
	}
	
	pseudo = c();
	for (trial in 1:trials)
	{
		pseudoData = sample(data, length(data), replace = T);
		pseudo = c(pseudo, eval(call(Func, pseudoData)));
	}
	
	half = (1 - interval) / 2;
	int = c(half, interval + half);
	conf.int = quantile(pseudo, int, na.rm=T);
	
	value = eval(call(Func, data[!is.na(data)]));
	if (plots)
	{
		par(mfrow = c(1, 2));
		hist(pseudo, 
			 breaks = length(pseudo)/100, 
			 col = col, border = border, 
			 main = paste('bootstrapped ', Func, '\n(n=', trials, ')', sep = ''), 
			 xlab = ''
			 );
		abline(v = value, col = line.col);
		abline(v = conf.int[1], col = line.col, lty = 'dashed');
		abline(v = conf.int[2], col = line.col, lty = 'dashed');
		stripchart(data, vertical = T, 
				   pch = pch, bg = col, 
				   frame.plot = F, 
				   main = paste(interval*100, '% conf.int around ', round(value, 1), 
				   				'\n[', round(conf.int[1], 1), ', ', 
				   				round(conf.int[2], 1), ']', 
				   				sep = ''
				   				)
				   );
		segments(.8, value, 1.2, value, col = line.col);
		text(1.26, value, Func);
		segments(.9, conf.int[1], 1.1, conf.int[1], col = line.col, lty = 'dashed');
		segments(.9, conf.int[2], 1.1, conf.int[2], col = line.col, lty = 'dashed');
	}
	
	
	output = list(value = value, conf.int = conf.int, dist = pseudo, parameters = match.call());
	
	return(output);
}





# combine datasets and check that data is numeric
# returns list of $group1, $group2
.validateBootstrapData = function(group1, group2) {
	
	if (!is.numeric(c(group1,group2))) stop('DATA CONTAINS NON-NUMERIC VALUES.');
	if (sum(is.na(c(group1,group2))) > 0)
	{
		group1 = group1[!is.na(group1)];
		group2 = group2[!is.na(group2)];
		warning('NAs removed from one or both groups, check your data');
	}
	return(list(group1=group1, group2=group2));
}


################################################################################################
#### bootstrap2independent() tests two numeric data vectors representing independent groups ####
################################################################################################

bootstrap2independent = function(x, y, 
						  		 trials = 10000, 
						  		 Func = 'mean', 
						  		 replace = T, 
						  		 plots = F, # TODO change back to T
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, #behaviorname
						  		 groupNames = c('group1', 'group2'), #control, CRISPR
						  		 jitter = .15, 
						  		 pch = 21,
						  		 boxLineMedian = F,
						  		 printResults = F,
						  		 verbose = F,
						  		 outfile = "__June/junkruns/0148" # NULL
						  		 )  #save as JPG: look at fishstudies/_code/bootstrapFunctions_6-16-13.R
{
	validateResult <- .validateBootstrapData(x, y);
	group1 = validateResult$group1;
	group2 = validateResult$group2;
	boxNULL = c(group1, group2);
	
	if (sum(boxNULL) == 0) {
		# BUG: If both groups are all 0, it complains. Could add a check here to fix that.
		# TODO basically return NA.
		# return(list(stat = eval(call(Func, group1)) - eval(call(Func, group2)), 
				# stat.reflect = reflect,
				# p.value = p, 
				# p.left = p_left,
				# p.right = p_right,
				# null.dist = statsNULL,
				# null.dist.mean = midNULL, 
				# replacement = replace,
				# data = data,
				# parameters = match.call()
				# );)
	}
	
	
	data = list(group1, group2);
	names(data) = c(groupNames[1], groupNames[2]);
	
	# compute statistic
	stat = eval(call(Func, group1)) - eval(call(Func, group2));
	
	# resample to create null distribution of statistic
	if (verbose)
	{
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: difference of group ', Func, 's\n', sep = '');
		cat(' Resampling with replacement: ', replace, '\n\n', sep = '');
	}
	statsNULL = c();
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trials <= 20000)
			{
				if (trial %% 1000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
			else if (trials > 20000)
			{
				if (trial %% 10000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
		}
		pseudo_group1 = sample(boxNULL, length(group1), replace = replace);
		pseudo_group2 = sample(boxNULL, length(group2), replace = replace);
		pseudo_stat_null = eval(call(Func, pseudo_group1)) - eval(call(Func, pseudo_group2));
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	
	# compute p-values
	midNULL = mean(statsNULL);
	reflect = midNULL - (stat - midNULL);
	if (stat < 0)
	{ 
		p_left = sum(statsNULL < stat) / trials;
		p_right = sum(statsNULL > reflect) / trials;
		p = p_left + p_right;
	} 
	else if (stat > 0)
	{
		p_right = sum(statsNULL > stat) / trials;
		p_left = sum(statsNULL < reflect) / trials;
		p = p_right + p_left;
	}
	else
	{
		stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
	}
	
	# histograms of datasets and null distribution, and boxplot with p-value
	if (plots)
	{
		# add n to group names
		if (!is.null(outfile)) jpeg(filename = outfile, width = 10, height = 10, units = "in", quality = 100, res = 150, type = "quartz");
		groupNames[1] = paste(groupNames[1], ' (n=', length(group1), ')', sep = '');
		groupNames[2] = paste(groupNames[2], ' (n=', length(group2), ')', sep = '');
		
		par(mfrow = c(2, 2));
		
		# compute histograms for each group
		hist_group1 = hist(group1, breaks = length(group1), plot = F);
		hist_group2 = hist(group2, breaks = length(group2), plot = F);
		
		# compute axis limits
		ymax = max(c(hist_group1$counts, hist_group2$counts));
		xmin = min(c(hist_group1$breaks, hist_group2$breaks));
		xmax = max(c(hist_group1$breaks, hist_group2$breaks));
		
		# check for dataDescriptor 
		if (!is.null(dataDescriptor) & is.character(dataDescriptor))
		{
			data.lab = dataDescriptor;
		}
		else
		{
			data.lab = '';
			warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
		}
		
		# plot group histograms with lines at means/medians
		plot(hist_group1, 
			 main = groupNames[1], 
			 ylim = c(0, ymax), xlim = c(xmin, xmax), 
			 col = col, border = border, 
			 xlab = data.lab
			 );
		abline(v = eval(call(Func, group1)), col = col.line);
		plot(hist_group2, 
			 main = groupNames[2], 
			 ylim = c(0, ymax), xlim = c(xmin, xmax), 
			 col = col, border = border, 
			 xlab = data.lab
			 );
		abline(v = eval(call(Func, group2)), col = col.line);
		
		# plot histogram of null distribution with lines at stat and reflected stat
		hist(statsNULL, 
		 	 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 	 xlab = paste('statistic (', Func, 's, replace = ', substr(replace, 1, 1), ')', sep = ''),
		 	 col = col, border = border,
		 	 breaks = floor(trials / 100)
		 	 );
		abline(v = stat, col = col.line);
		abline(v = reflect, col = col.line, lty = 'dashed');
		
		# boxplot of groups with individual data points
		toPlot = boxNULL;
		grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], length(group2)));
		
		# check if line in box should be mean or median (default)
		if (boxLineMedian)
		{
			medlty = 'solid';
		}
		else
		{
			medlty = 'blank';
		}
		
		if (p == 0)
		{
			toPaste = paste('p < 1e-5', sep = '');
		}
		else
		{
			toPaste = paste('p = ', signif(p, 2), sep = '');
		}
		
		# draw boxes
		boxplot(toPlot ~ grp,
				ylab = data.lab,
				medlty = medlty,
				main = toPaste
				);
		
		# draw mean lines if specified
		if (!boxLineMedian)
		{
			segments(0.6, mean(group1), 1.4, mean(group1), lwd = 2);
			segments(1.6, mean(group2), 2.4, mean(group2), lwd = 2);
		}
		
		# add data points
		stripchart(toPlot ~ grp,
				   vertical = T,
				   add = T,
				   method = 'jitter',
				   jitter = jitter,
				   pch = pch,
				   bg = col,
				   cex = 1.5
				   );
		if (!is.null(outfile)) dev.off();
	}
	
	# build output list
	output = list(stat = stat, 
				stat.reflect = reflect,
				p.value = p, 
				p.left = p_left,
				p.right = p_right,
				null.dist = statsNULL,
				null.dist.mean = midNULL, 
				replacement = replace,
				data = data,
				parameters = match.call()
				);
	
	# print summary of results to screen if specified			
	if (printResults)
	{
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist and $parameters\n');
	}
	
	return(output);
}

################################################
#### bootstrap2paired  ####
################################################

bootstrap2paired = function(condition1, condition2, 
						  		 trials = 10000, 
						  		 Func = 'mean',
						  		 plots = T, 
						  		 abs.diffs = F,
						  		 col = 'grey',
						  		 border = 'darkgrey',
						  		 col.line = 'red',
						  		 dataDescriptor = NULL, 
						  		 conditionNames = c('condition1', 'condition2'), 
						  		 pch = 21,
						  		 cex.lab = 1.2, cex.axis = 1.2,
						  		 printResults = T,
						  		 verbose = T,
						  		 ...
						  		 )
{
	# check that data is numeric
	numCheck1 = is.numeric(condition1);
	numCheck2 = is.numeric(condition2);
	if (!numCheck1 | !numCheck2)
	{
		stop('DATA CONTAINS NON-NUMERIC VALUES...\n   I ONLY EAT NUMBERS!!!\n    GIVE ME NUMBERS!!!!!');
	}
	
	# check that conditions have same number of data points
	if (length(condition1) != length(condition2))
	{
		stop('GROUPS ARE DIFFERENT SIZES...\n   THIS DATA SHOULD BE PAIRED!!!');
	}
	
	# check for missing data
	NAcheck1 = is.na(condition1);
	NAcheck2 = is.na(condition2);
	removeMe = NAcheck1 | NAcheck2;
	if (sum(removeMe) > 0)
	{
		condition1 = condition1[!removeMe];
		condition2 = condition2[!removeMe];
		warning('NAs in one/both conditions, corresponding data points removed from BOTH');
	}
	
	# compute statistic
	diffs = condition1 - condition2; 
	stat = eval(call(Func, diffs));
	if (abs.diffs)
	{
		stat = abs(stat);
	}
	
	# store data for output
	data = data.frame(condition1, condition2, diffs);
	names(data) = c(conditionNames[1], conditionNames[2], 'diffs');
	
	# resample to create null distribution of statistic
	if (verbose)
	{
		cat('...........................................\n');
		cat('Testing: ', dataDescriptor, '\n\n', sep = '');
		cat('Building null distribution...\n');
		cat(' Test statistic: ', Func, ' of individual differences\n\n', sep = '');
	}
	statsNULL = c();
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trials <= 20000)
			{
				if (trial %% 1000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
			else if (trials > 20000)
			{
				if (trial %% 10000 == 0)
				{
					cat('  Run ', trial, '\n', sep = '');
				}
			}
		}
		signs = sample(c(-1, 1), length(diffs), replace = T);
		pseudo_diffs = signs * diffs;
		pseudo_stat_null = eval(call(Func, pseudo_diffs));
		if (abs.diffs)
		{
			pseudo_stat_null = abs(pseudo_stat_null);
		}
		statsNULL = c(statsNULL, pseudo_stat_null);
	}
	
	# compute p-values
	if (abs.diffs)
	{
		p = sum(statsNULL > stat) / trials;
	}
	else if (!abs.diffs)
	{
		midNULL = mean(statsNULL);
		reflect = midNULL - (stat - midNULL);
		if (stat < 0)
		{ 
			p_left = sum(statsNULL < stat) / trials;
			p_right = sum(statsNULL > reflect) / trials;
			p = p_left + p_right;
		} 
		else if (stat > 0)
		{
			p_right = sum(statsNULL > stat) / trials;
			p_left = sum(statsNULL < reflect) / trials;
			p = p_right + p_left;
		}
		else
		{
			stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
		}
	}
	
	if (plots)
	{
		# check for dataDescriptor 
		if (!is.null(dataDescriptor) & is.character(dataDescriptor))
		{
			data.lab = dataDescriptor;
		}
		else
		{
			data.lab = '';
			warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
		}
		
		par(mfrow = c(1, 3));
		
		# compute histogram of individual differences
		hist.diffs = hist(diffs, breaks = length(diffs), plot = F);
		xmax = max(abs(hist.diffs$breaks));
	
		# plot histogram
		plot(hist.diffs,
			 xlim = c(-xmax, xmax),
			 main = paste('Individual differences with ', Func, '\n(n=', length(diffs), ')', sep = ''),
			 xlab = paste(dataDescriptor, ' (', conditionNames[1], ' - ', conditionNames[2], ')', sep = ''),
			 col = col, border = border,
			 cex.lab = cex.lab,
			 cex.axis = cex.axis
			 );
		abline(v = mean(diffs), col = col.line);
		
		# plot histogram of null distribution with lines at stat and reflected stat
		if (abs.diffs)
		{
			xlabNULL = paste('statistic: abs(', Func, ' of individual differences)', sep = '');
		}
		else if (!abs.diffs)
		{
			xlabNULL = paste('statistic: ', Func, ' of individual differences', sep = '');
		}
		hist(statsNULL, 
		 	 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''), 
		 	 xlab = xlabNULL,
		 	 col = col, border = border,
		 	 breaks = floor(trials / 100),
		 	 cex.lab = cex.lab,
			 cex.axis = cex.axis
		 	 );
		if (abs.diffs)
		{
			abline(v = stat, col = col.line);
		}
		else if (!abs.diffs)
		{
			abline(v = stat, col = col.line);
			abline(v = reflect, col = col.line, lty = 'dashed');
		}
	
		# boxplot of conditions with individual data points
		toPlot = data[, 1:2];
		
		if (p == 0)
		{
			toPaste = paste('p < 1e-5 (n=', length(diffs), ')', sep = '');
		}
		else
		{
			toPaste = paste('p = ', signif(p, 2), ' (n=', length(diffs), ')', sep = '');
		}
		
		# draw boxes
		boxplot(toPlot,
				ylab = data.lab,
				#medlty = medlty,
				main = toPaste,
				frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,
				...
				);
				
		# add data points
		stripchart(toPlot,
				   vertical = T,
				   add = T,
				   pch = pch,
				   bg = col,
				   cex = 1.5
				   );
				   
		# connect paired points across conditions
		for (row in 1:nrow(data))
		{
			segments(1, data[row, 1], 2, data[row, 2], col = border);
		}
		
	}
	
	# build output list
	if (abs.diffs)
	{
		output = list(stat = stat,
					  p = p,
					  null.dist = statsNULL,
					  data = data,
					  parameters = match.call()
					  );
	}
	else if (!abs.diffs)
	{
		output = list(stat = stat, 
					  stat.reflect = reflect,
					  p = p, 
					  p.left = p_left,
				  	  p.right = p_right,
				  	  null.dist = statsNULL,
				  	  null.dist.mean = midNULL,
				  	  data = data,
				  	  parameters = match.call()
				  	  );
	}
				  
	# print summary of results to screen if specified			
	if (printResults)
	{
		cat('\n');
		print(lapply(output, head));
		cat('   ... only showing first few values of $null.dist, $data, and $parameters\n');
	}

	return(output);
}

################################################################################################
####  ####
################################################################################################

calculateF = function(data = 'matrix or dataframe with groups in columns', 
					  Func = 'mean',
					  absDiffs = T
					  )
{
	# check data
	if (!is.data.frame(data) & !is.matrix(data))
	{
		stop('Data is not stored in matrix or dataframe...\n   FIX IT!');
	}
	if (mode(data) != 'numeric')
	{
		stop('DATA ISN\'T NUMERIC...\n   I ONLY EAT NUMBERS!!!');
	}
	if (nrow(data) < ncol(data))
	{
		warning('More groups than data points per group...\n   ARE YOU SURE?');
	}
	data = as.data.frame(data);
	
	# group info
	nGroups = ncol(data);
	groupNums = apply(data, 2, get(Func), na.rm = T);
	groupNAs = apply(apply(data, 2, is.na), 2, sum);
	groupLengths = nrow(data) - groupNAs;	
	# store for output
	groupInfo = list(data = data, 
					 nGroups = nGroups, 
					 groupNums = groupNums, 
					 groupNAs = groupNAs, 
					 groupLengths = groupLengths
					 );		
	
	# grand mean/median/Func
	grandTop = sum(groupLengths * groupNums);
	grandBot = sum(groupLengths);
	grandNum = grandTop / grandBot;
	# store for output
	grandNumInfo = list(grandTop = grandTop, 
						grandBot = grandBot, 
						grandNum = grandNum
						);
	
	# test statistic numerator
	diffs = grandNum - groupNums;
	if (absDiffs)
	{
		diffs = abs(diffs);
	}
	numerator = sum(groupLengths * diffs);
	# store for output
	Fnum = list(diffs = diffs,
				numerator = numerator
				);
	
	# test statistic denominator
	toSum = list();
	sums = c();
	for (gp in 1:nGroups)
	{
		temp = data[, gp] - groupNums[gp];
		temp = temp[!is.na(temp)];
		if (absDiffs)
		{
			temp = abs(temp);
		}
		toSum[[gp]] = temp;
		sums[gp] = sum(toSum[[gp]]);
	}
	denominator = sum(sums);
	# store for output
	Fden = list(toSum = toSum,
				sums = sums,
				denominator = denominator
				);
	
	Fstat = numerator / denominator; 
	# store for output
	FstatInfo = list(Fnum = Fnum,
					 Fden = Fden,
					 Fstat = Fstat
					 );
	
	out = list(groupInfo = groupInfo, 
			   grandNumInfo = grandNumInfo,
			   FstatInfo = FstatInfo,
			   parameters = match.call()
			   );
	
	return(out);
}

#######################################
####  ####
#######################################

bootstrapANOVA = function(data = 'matrix or dataframe with groups in columns', 
					  	  Func = 'mean',
					  	  absDiffs = T,
					  	  trials = 10000, 
					  	  replace = F,
					  	  plots = T,
					  	  groupNames = NULL,
					  	  pch = 21,
					  	  col = 'grey', border = 'darkgrey',
					  	  line.col = 'red', 
					  	  verbose = T
					  	  )
{
	# calculate actual statistic and get group info
	calcF = calculateF(data = data, Func = Func, absDiffs = absDiffs);
	Fstat = calcF$FstatInfo$Fstat;
	groupLengths = calcF$groupInfo$groupLengths;#print(groupLengths)
	nGroups = calcF$groupInfo$nGroups;
	
	# combine groups
	boxNULL = unlist(calcF$groupInfo$data);
	boxNULL = boxNULL[!(is.na(boxNULL))];#print(boxNULL)
	statsNULL = c();
	
	for (trial in 1:trials)
	{
		if (verbose)
		{
			if (trial %% 1000 == 0)
			{
				cat('  Run ', trial, '\n', sep = '');
			}
			
		}
		pseudoNULL = matrix(nrow = max(groupLengths), ncol = nGroups);		
		for (gp in 1:ncol(pseudoNULL))
		{
			if (groupLengths[gp] < nrow(pseudoNULL))
			{
				numNA = nrow(pseudoNULL) - groupLengths[gp];
				pseudoNULL[1:numNA, gp] = NA;
				pseudoNULL[(numNA+1):nrow(pseudoNULL), gp] = sample(boxNULL, groupLengths[gp], replace = replace);
			}
			else
			{
				pseudoNULL[, gp] = sample(boxNULL, groupLengths[gp], replace = replace);
			}	
		}
		#print(pseudoNULL)
		pseudoF = calculateF(data = pseudoNULL, Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat;
		statsNULL = c(statsNULL, pseudoF);
	}
	
	# compute p-values
	if (absDiffs)
	{
		p = sum(statsNULL > Fstat) / trials;
		p_left = NULL;
		p_right = NULL;
	}
	else
	{
		midNULL = mean(statsNULL);
		reflect = midNULL - (Fstat - midNULL);
			if (Fstat < 0)
			{
				p_left = sum(statsNULL < Fstat) / trials;
				p_right = sum(statsNULL > reflect) / trials;
				p = p_left + p_right;
			}
			else if (Fstat > 0)
			{
				p_right = sum(statsNULL > Fstat) / trials;
				p_left = sum(statsNULL < reflect) / trials;
				p = p_right + p_left;
			}
			else
			{
				stop('Either statistic==0 or something else is wrong...\n  Figure it out or find Austin!');
			}
	}
	

	
	if (plots)
	{
		
		if (verbose)
		{
			cat('  Computing 95% confidence intervals\n');
		}
		values = c();
		ints = matrix(nrow = ncol(data), ncol = 2);
		for (gp in 1:ncol(data))
		{
			temp = confidenceInterval(data[, gp], plots = F, Func = Func);
			values = c(values, temp$value);
			ints[gp, ] = temp$conf.int;
		}
		print(values)
		par(mfrow = c(1, 2));
		hist(statsNULL, 
			 xlim = c(min(statsNULL), max(statsNULL, Fstat)), 
			 main = paste('bootstrap null distribution\n(n=', trials, ')', sep = ''),
			 col = col, border = border, 
			 xlab = ''
			 );
		if (absDiffs)
		{
			abline(v = Fstat, col = line.col);
		}
		else 
		{
			abline(v = Fstat, col = line.col);
			abline(v = reflect, col = line.col, lty = 'dashed');
		}
		
		if (is.null(groupNames))
		{
			groupNames = colnames(data);
		}
		else if (!is.null(groupNames) & length(groupNames)==ncol(data))
		{
			groupNames = groupNames;
		}
		boxplot(data, 
			    names = groupNames, 
			    main = paste('p = ', p, sep = ''), 
			    border = 'lightgrey', 
			    medlty = 'blank', 
			    boxwex = .5
			    );
			
		stripchart(as.data.frame(data), vertical = T, add = T, pch = pch, bg = col);
		
		coMat = matrix(nrow = ncol(data), ncol = 2);
		coMat[1, ] = c(.75, 1.25);
		for (gp in 2:nrow(coMat))
		{
			coMat[gp, ] = coMat[gp - 1, ] + 1;
		}
		for (gp in 1:length(values))
		{
			segments(coMat[gp, 1], values[gp], 
				     coMat[gp, 2], values[gp], 
				     lwd = 1, col = line.col
				     );
			segments(coMat[gp, 1], ints[gp, 1], 
					 coMat[gp, 2], ints[gp, 1], 
					 lwd = 1, lty = 'dashed', col = line.col
					 );
			segments(coMat[gp, 1], ints[gp, 2], 
				     coMat[gp, 2], ints[gp, 2], 
				     lwd = 1, lty = 'dashed', col = line.col
				     );
		}
		
	}
	
	output = list(stat = Fstat,
			      data = data,
			      null.dist = statsNULL,
			      p = p,
			      p.left = p_left,
			      p.right = p_right,
			      calcF = calcF,
			      conf.ints = ints,
			      parameters = match.call()
			      );
	
	return(output);
}
			