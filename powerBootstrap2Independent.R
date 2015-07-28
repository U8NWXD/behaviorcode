############################################################################################
### function to compute power of a bootstrap comparison between two independent datasets ###
############################################################################################
### power is an estimate of how often a significant difference would be detected if it actually exists, usually want power > 0.8
### it depends on group size, signifiance threshold, and effect size
###   here, group size is inferred from input data, could add option to resample input data for purpose of simulating larger groups
###         significance threshold is hard coded as 95th%tile, could add arg for user to test different cutoffs
###         effect size is computed from input data if both 'ctrl' and 'exp' data are provided
###         effect size is defined by user with 'change' if only 'ctrl' is provided
###           e.g. change=0.8 and change=1.5 will simulate a 20% decrease and a 50% increase from ctrl, repsectively
###           could add options for more user control over characteristics of simulated group, e.g. size
###           
### if both 'ctrl' and 'exp' are provided will compute retrospective power
### if only 'ctrl' provided function tries to simulate an 'exp' group and compute prospective power for a certain effect size (defined by 'change')
###   if 'exp' is NULL and 'change' is not numeric function will throw error       

### 'ctrl' and 'exp' (if provided) are numeric vectors, don't need to be equal lengths
### 'change' (if provided) is a double
### 'Func' is string naming function to compute a descriptor of each dataset, e.g. mean, median, etc
###   test statistic will be computed as the difference between the descriptor values for ctrl and exp groups
### 'trials' is integer to set the number of resampling runs
### 'abs' is boolean indicating whether to ignore the sign of the test statistic
###   in practice this usually makes little (if any) difference to the power calculation
### 'plot' is boolean indicating whether to plot histograms of output distributions
### 'breaks' (if provided) is integer defining number of bins for histograms
### 'cex.main' sets sizes of plot and histogram titles
### 'verbose' is boolean indicating whether to print trial number every 1000 trials 

powerBootstrap2Independent = function (ctrl, exp=NULL, change=NULL, Func='mean', trials=10000, abs=F, plot=T, breaks=NULL, cex.main=1.1, verbose=T, outfile = NULL, sigThreshold = .05, ...) {
	# check if only ctrl group is provided
	if (is.null(exp)) {
		# if yes, simulate exp group
		if (is.numeric(change)) {
			exp = change * ctrl;
		# if yes, but no value is provided for 'change' throw error
		} else {
			stop('if no "exp" data provided "change" must be numeric');
		}
	}
	
	# remove NAs
	if (sum(is.na(c(ctrl,exp))) > 0) {
		ctrl = ctrl[!is.na(ctrl)];
		exp = exp[!is.na(exp)];
		warning('NAs removed from one or both groups, check your data. Power test output will be based on number of not-NA data points, NOT the group size provided.', immediate. = T);
	}
	
	# compute actual value of statistic (stat) and %change (delta) from ctrl to exp group 
	old = eval(call(Func, ctrl));
	new = eval(call(Func, exp));
	stat = old - new;
	if (abs) { stat = abs(stat) }
	delta = (new - old) / old;
	
	# combine ctrl and exp for resampling to create null distribution
	boxNULL = c(ctrl, exp);
	# initialize vecs to hold values of stat for null and alt distributions 
	statsNULL = vector(mode='numeric',length=trials);
	statsALT = vector(mode='numeric',length=trials);
	for (trial in 1:trials) {
		if (verbose & (trial %% 1000 == 0)) { cat(paste(' ',trial,sep='')) } 
		
		# compute stat for null distribution by resampling from combined datasets 
		pseudo_ctrl = sample(boxNULL, length(ctrl), replace=T);
		pseudo_exp = sample(boxNULL, length(exp), replace=T);
		pseudo_stat_null = eval(call(Func, pseudo_ctrl)) - eval(call(Func, pseudo_exp));
		statsNULL[trial] = pseudo_stat_null;
		
		# compute stat for alt distribution by resampling within each dataset
		pseudo_ctrlALT = sample(ctrl, length(ctrl), replace=T);
		pseudo_expALT = sample(exp, length(exp), replace=T);
		pseudo_stat_alt = eval(call(Func, pseudo_ctrlALT)) - eval(call(Func, pseudo_expALT));
		statsALT[trial] = pseudo_stat_alt;
	}
	
	# get power by computing how much of alt distribution is more extreme than 95th%tile of null distribution
	if (abs) {
		statsNULL = abs(statsNULL);
		NULLconfint = quantile(statsNULL, 1 - sigThreshold, na.rm=T); # TODO not a magic number
		statsALT = abs(statsALT);
		right_count = sum(statsALT > NULLconfint,na.rm=T);
		left_count = NULL;
		power = right_count / trials;
	} else {
		NULLconfint = quantile(statsNULL, c(sigThreshold/2, 1 - sigThreshold / 2), na.rm=T);
		right_count = sum(statsALT > NULLconfint[2],na.rm=T);
		left_count = sum(statsALT < NULLconfint[1],na.rm=T);
		power = (right_count + left_count) / trials;
	}
	
	# plot histograms of alt and null distributions with vertical line(s) at significance cutoff(s)
	if (plot) {
		if (!is.null(outfile)) jpeg(filename = outfile, width = 6, height = 10, units = "in", quality = 100, res = 150, type = "quartz");
		if (is.null(breaks)) { breaks = trials / 2 }
		par(mfrow=c(2, 1), oma=c(0,0,3,0), mar=c(2,4,3,2));
		# print(breaks)
		# print(head(statsNULL))
		histNULL = hist(statsNULL, breaks=breaks, plot=F);
		histALT = hist(statsALT, breaks=breaks, plot=F);
		
		# get extreme valuess of axes to line up alt and null histograms
		ymax = max(c(histNULL$counts, histALT$counts));
		xmax = max(c(histNULL$breaks, histALT$breaks));
		xmin = min(c(histNULL$breaks, histALT$breaks));
		
		plot(histNULL, main='null dist', col='grey', border='darkgrey', xlab='', ylim=c(0, ymax), xlim=c(xmin, xmax), cex.main=cex.main*.85, ...);
		abline(v = NULLconfint, col = 'red', lty = 'dashed');
		plot(histALT, main='alternate dist', col='grey', border='darkgrey', xlab='', ylim=c(0, ymax), xlim=c(xmin, xmax), cex.main=cex.main*.85, ...);
		abline(v = NULLconfint, col = 'red', lty = 'dashed');
		title(main=paste('Test statistic: difference of group ', Func, 's\npower = ', signif(power, 2), sep=''), outer=T, cex.main=cex.main);
		if (!is.null(outfile)) dev.off()
	}
	
	return(list(ctrl=ctrl, exp=exp, stat=stat, delta=delta, confint=NULLconfint, rightcount=right_count, leftcount=left_count, power=power, null.dist=statsNULL, alt.dist=statsALT));
}