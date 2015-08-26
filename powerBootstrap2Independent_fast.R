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
### 'change' (if provided) is a double (default=NULL)
### 'Func' is string naming function to compute a descriptor of each dataset, e.g. mean, median, etc (default='mean')
###   test statistic will be computed as the difference between the descriptor values for ctrl and exp groups
### 'trials' is integer to set the number of resampling runs (default=10000)
### 'abs' is boolean indicating whether to ignore the sign of the test statistic (default=FALSE)
###   in practice this usually makes little (if any) difference to the power calculation
### 'plot' is boolean indicating whether to plot histograms of output distributions, jpeg will be saved to working directory if plotFile=NULL (default=TRUE)
### 'breaks' (if provided) is integer defining number of bins for histograms (default=NULL)
### 'cex.main' sets sizes of plot and histogram titles (default=1.1)
### 'verbose' is boolean indicating whether to print trial number every 1000 trials (default=TRUE)
### 'plotFile' is string for providing a filename base for jpeg of results (if plot=TRUE), will be appended by '_power.jpeg' (default=NULL)
###   if NULL, plotFile will be set to the name of the 'ctrl' group and saved in the working directory 
### 'increase_n' is integer defining a number of additional subjects to simulate in each group via extra resampling (default=0)
###   e.g. if increase_n=2 and the actual group sizes are ctrl=5 and exp=3, pseudo-group sizes during the power analysis will be ctrl=7 and exp=5

powerBootstrap2Independent = function (ctrl, exp=NULL, change=NULL, Func='mean', trials=10000, abs=FALSE, plot=TRUE, breaks=NULL, cex.main=1.1, verbose=TRUE, plotFile=NULL, increase_n=0, ...) {
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
	
	if (sum(is.na(c(ctrl, exp)))) {
		warning("NAs present in data. Removing them and adjusting group sizes down as a result. You may want to set the increase_n parameter to the number of NAs there were to preserve group sizes.")
		ctrl = ctrl[!is.na(ctrl)]
		exp = exp[!is.na(exp)]
	}
	
	# compute actual value of statistic (stat) and %change (delta) from ctrl to exp group 
	old = eval(call(Func, ctrl, na.rm=T));
	new = eval(call(Func, exp, na.rm=T));
	stat = old - new;#print(stat)
	if (stat == 0 | length(ctrl[!is.na(ctrl)]) < 3 | length(exp[!is.na(exp)]) < 3) {
		warning(paste('Either no difference between group ', Func, 's or one of the groups has < 3 valid data points, thus power is undefined. \nReturning a list with one element: $power=NA', sep=''));
		return(list(power=NA));
	}
	if (abs) { stat = abs(stat) }
	delta = (new - old) / old;
	
	# combine ctrl and exp for resampling to create null distribution
	boxNULL = c(ctrl, exp);
	# initialize vecs to hold values of stat for null and alt distributions 
	statsNULL = vector(mode='numeric',length=trials);
	statsALT = vector(mode='numeric',length=trials);
	fxn = get(Func)
	ctrlsamplelen = length(ctrl)+increase_n;
	expsamplelen = length(exp)+increase_n;
	for (trial in 1:trials) {
		if ((trial %% 1000 == 0) & verbose) { cat(paste(' ',trial,'\n',sep='')) } 
		
		# compute stat for null distribution by resampling from combined datasets 
		# pseudo_ctrl = sample(boxNULL, ctrlsamplelen, replace=T);
		# pseudo_exp = sample(boxNULL, expsamplelen, replace=T);
		# pseudo_stat_null = fxn(pseudo_ctrl, na.rm=T) - fxn(pseudo_exp, na.rm=T);
		# statsNULL[trial] = pseudo_stat_null;
		statsNULL[trial] = fxn(sample(boxNULL, ctrlsamplelen, replace=T)) - fxn(sample(boxNULL, expsamplelen, replace=T));
		
		# compute stat for alt distribution by resampling within each dataset
		# pseudo_ctrlALT = sample(ctrl, ctrlsamplelen, replace=T);
		# pseudo_expALT = sample(exp, expsamplelen, replace=T);
		# pseudo_stat_alt = fxn(pseudo_ctrlALT, na.rm=T) - fxn(pseudo_expALT, na.rm=T);
		# statsALT[trial] = pseudo_stat_alt;
		statsALT[trial] = fxn(sample(ctrl, ctrlsamplelen, replace=T)) - fxn(sample(exp, expsamplelen, replace=T));
	}
	
	# get power by computing how much of alt distribution is more extreme than 95th%tile of null distribution
	if (abs) {
		statsNULL = abs(statsNULL);
		NULLconfint = quantile(statsNULL, 0.95, na.rm=T);
		statsALT = abs(statsALT);
		right_count = sum(statsALT > NULLconfint,na.rm=T);
		left_count = NULL;
		power = right_count / trials;
	} else {
		NULLconfint = quantile(statsNULL, c(0.025, 0.975), na.rm=T);
		right_count = sum(statsALT > NULLconfint[2],na.rm=T);
		left_count = sum(statsALT < NULLconfint[1],na.rm=T);
		power = (right_count + left_count) / trials;
	} # TODO get opts I added to the old power test and put them hereee.
	
	# plot histograms of alt and null distributions with vertical line(s) at significance cutoff(s)
	if (plot) {
		if (is.null(plotFile)) { 
			plotFile = paste(deparse(substitute(ctrl)), '_power.jpg', sep='');
		} else if (!is.character(plotFile)) { 
			stop('Arg plotFile must be NULL or character');
		}
		jpeg(file=plotFile, width=8, height=8, units='in', quality=100, type='quartz', res=150);
		if (is.null(breaks)) { breaks = trials / 2 }
		par(mfrow=c(2, 1), oma=c(0,0,3,0), mar=c(2,4,3,2));
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
		title(main=paste('Test statistic: difference of group ', Func, 's, ', trials, ' runs\npower = ', signif(power, 2), sep=''), outer=T, cex.main=cex.main);
		dev.off();
	}
	
	return(list(ctrl=ctrl, exp=exp, stat=stat, delta=delta, confint=NULLconfint, rightcount=right_count, leftcount=left_count, power=power, null.dist=statsNULL, alt.dist=statsALT));
}