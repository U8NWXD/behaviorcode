library(stringr)
options(stringsAsFactors = FALSE)
source("~/Desktop/Katrina/behavior_code/bootstrap_tests_June2013_STABLE.R")
#use color=blue in .dot output script to make separate sets of lines for 1st follower, 2nd, etc.
# Find a way to represent AB -> C, ABC -> D instead of only A -> B probabilities
# collapse statistics across animals. Arrows only allowed from subj to different subj.
# frame counts!!!
# Ask scott - what defines start of spawning???
# maybe paths leading up to spawning - size of circle reps # of paths. Forks and such.
#heatmaps with k!
# HEATMAP OF COUNTS FOR INDIVIDUAL SAMPLES. 6x6. DO THE ONES CLUSTER.
#bouts
# TODO paired observations instead of two groups of fish

# behaviors <- names(table(sjan_data[[6]]$dat$beh))
# behaviors <- c(behaviors[1:2], 'c', behaviors[3:length(behaviors)]) #hacky and not generally correct
# mat = matrix(nrow = length(behaviors), ncol = 6)
# rownames(mat)<- behaviors
# colnames(mat)<- names(sjan_data)
# for (i in 1:length(sjan_data)) {tab = table(sjan_data[[i]]$dat$beh);
# for (j in 1:length(behaviors)) { if (behaviors[j] %in% names(tab)) {mat[behaviors[j],i] <- tab[behaviors[j]]} else {mat[behaviors[j],i] <- 0}
# }}

# > eval(parse(text="tmp[order(tmp[,4], tmp[,3], tmp[,2], tmp[,6], tmp[,7]),]"))



# heatmap(cor(t(probma)), symm = TRUE)



# TODO add option for all compare fxns verbose = false

# TODO add or ask about folder slash-at-end. Maybe do this in the Big Shell that asks for ONE outfile path.
# TODO make the Big Shell

#####################################################################################################
## TINY HELPERS                                                                                    ##
#####################################################################################################

# Checks that <colName> is either a valid column index for <object> (a data frame or
# matrix) or a valid column name for <object>.
.checkColumnOK = function(colName, object) {
	return((is.numeric(colName) && colName %% 1 == 0 && colName > 0 && colName <= dim(object)[2]) ||
		   (is.character(colName) && colName %in% dimnames(object)[[2]]))
}

# Sorts df by column time and returns the result.
.sortByTime = function(df) {
	df = df[order(as.numeric(df$time)),];
	dimnames(df)[[1]] <- 1:(length(df$time));
	return(df);
}

# Converts a time in the format "MINUTES:SECONDS" to a numver of seconds.
# Fractional seconds are ok (e.g. "3:14.15")
.timeToSeconds = function(clockTime) {
	nMinutes = as.numeric(gsub(":.*$", "", clockTime));
	nSeconds = as.numeric(gsub("^.*:", "", clockTime));
	return(60 * nMinutes + nSeconds);
}

# Gets user input as either a yes or no answer to <prompt>.
# If the inputs starts with 'y' or 'Y', returns true; 'n' or 'N' returns false; other letters result in a reprompt.
.getYesOrNo = function(prompt) {
	response = tolower(substr(readline(prompt), 1, 1));
	while (response != 'y' && response != 'n') response = tolower(substr(readline("Please enter yes (y) or no (n): "), 1, 1));
	return(response == 'y');
}

# "Autocompletes" <input> to be one of the strings in <choices> by looking to see if there
# is exactly one string in <choices> that has <input> as a prefix.
# If case sensitivity is turned off (default), <input> will be modified to match the case
# of its match in <choices>.
.autocomplete = function(input, choices, caseSensitive = FALSE) {
	greplResult = if (caseSensitive) grepl(paste("^", input, sep = ""), choices) else grepl(paste("^", tolower(input), sep = ""), tolower(choices));
	if (sum(greplResult) == 1) {return(choices[greplResult]);}
	else if (!caseSensitive && tolower(input) %in% tolower(choices)) {return(choices[which(tolower(choices) == tolower(input))]);}
	else return(input);
}

#####################################################################################################
## READING DATA FROM SCORE LOGS                                                                    ##
#####################################################################################################

# Source: ethograms_from_scorevideo.R
# Calls .getData() on every file in <folderPath> and returns the results in a list.
# If groups = TRUE, will recurse through directories. Use with one folder for control animals,
# one folder for experimental condition 1, etc.
# Advanced: Assay start can be FALSE to turn off assay starts, or NULL or a vector of default assay starts
.getDataBatch = function (folderPath, groups = FALSE, assayStart = NA) {
	filenames = if (groups) {paste(folderPath,list.files(folderPath, pattern = "txt$", recursive = TRUE),sep = "");}
	            else {paste(folderPath,list.files(folderPath, pattern = "txt$"),sep = "");}
	data = list();
	if (!is.null(assayStart) && is.na(assayStart)) assayStart = if (.getYesOrNo("Did you mark assay starts in your score logs? ")) NULL else FALSE;
	for (f in 1:length(filenames)) {
		cat("Loading file \"", filenames[f], "\"...\n", sep = "");
		datOut = .getData(filenames[f], assayStart, single = FALSE);
		data[[f]] = datOut[[1]];
		assayStart = datOut[[2]];
		names(data)[f] = gsub(folderPath, "", filenames[f]);
	}
	
	# .printFindDupBehaviors(data);
	data <- .removeNALogs(data); # TODO is this what should be done about this?
	# print(lapply(data, function(f) {names(table(f$behavior))}));
	
	cat("Behaviors found:\n");
	.printFindDupBehaviors(data);
	cat('There may be some behaviors in the list above that should be combined (for example, "Female Follows" and "female follows").\n');
	if (.getYesOrNo("Are there any behaviors in the list that should be combined? ")) {
		data <- .promptToElimDups(data);
	}
	
	if(.getYesOrNo("Were all of your assays the same length of time? ")) {
		userInput = readline("Please enter the length of your assay in seconds.\n> ");
		while (!grepl("^[0-9]*$", userInput)) userInput = readline("Unreadable time format.\nPlease enter the length of your assay in seconds.\n> ");
		data = .filterDataList(data, endTime = as.numeric(userInput));
	}
	
	# TODO common errors check
	cat("Data entry complete.\n")
	return(.filterDataList(data, renameSubjects = TRUE));
}


# Reads in a score log (.txt file) and returns a data frame with columns
#   time    behavior    subject    type    pair_time    duration
# This data frame has an attribute called assay.start that gives the name
#   of the mark the user indicated was the assay start and the time in the
#   video that this mark occured. Regardless of the time given in this attribute,
#   the assay start is defined as time 0 in the timescale used for "time" and
#   "pair_time" in the data frame itself (that is, behavior times are given
#   in seconds after the assay start, NOT seconds after the video start).
# If called with single = FALSE (should only be done from .getDataBatch), this
#   is returned as the first element of a list whose second element is the updated
#   assayStart (since no pass-by-reference).
.getData = function (filename, assayStart = NULL, single = TRUE) { # TODO rewrite w/o single and with assayStart passed-by-ref in C.
	data0 = read.table(filename, fill=T, colClasses='character', sep='\t', header=F, quote='', blank.lines.skip=T, strip.white=T);
	desc_table = .getDescriptionTableFromRawData(data0);
	
	fpsRow = which(grepl('^Frames/sec of files:', data0[,1]));
	framesPerSecond = as.numeric(gsub('[^0-9]', '', data0[fpsRow,1]));
	
	asout = .getAssayStart(data0, assayStart, filename);
	startTime = asout$time;
	
	df = .parseFullLog(data0, desc_table, framesPerSecond);
	if (!is.data.frame(df)) {
		warning("EMPTY SCORE LOG.", immediate. = TRUE);
		return(if(single) df else list(df, asout$assayStart));
	} else if (length(df$behavior) < 10) {
		warning("Fewer than 10 behaviors in log.", immediate.=T)
	}
	
	
	df = .sortByTime(df);
	
	if (df$time[1] < startTime) warning("Some behavior(s) were scored before the assay start.", immediate. = TRUE);
	df$time <- df$time - startTime;
	df$pair_time <- df$pair_time - startTime;
	attr(df, 'assay.start') <- list(mark = asout$name, time = startTime);
	
	if (nrow(df)<1) {
		warning(paste('No data in ', filename, sep=''));
		return(list(data=df, codes=NULL));
	}
#	desc_table = desc_table[as.character(desc_table[, 1]) %in% as.character(df$behavior), ];

#   could return desc_table here if desired
	return(if (single) df else list(df, asout$assayStart));
}

# Helper function for .getData
# Gets the assay start mark with assistence from the user and returns a list whose first element
#   is the assay start time (or 0 if no assay start was marked) and whose second element is a
#   character vector of default assay starts.
.getAssayStart = function(data0, assayStart, logname) {
	if (is.logical(assayStart) && !assayStart) return(list(time = 0, name = NA, assayStart = assayStart));
	
	start_marks = which(data0=='MARKS') + 4;
	end_marks = length(data0[,1]) - 1;
	
	marks = strsplit(data0[start_marks:end_marks, ], '    ');
	marks = lapply(marks, function(f) gsub('^ *','', gsub(' *$','', f)));
	markNames = unlist(lapply(marks, function(f) f[[3]]));
	
	if (!is.null(assayStart) && sum(markNames %in% assayStart) == 1) {
		markIndex = which(markNames %in% assayStart);
		return(list(time = .timeToSeconds(marks[[markIndex]][2]), name = markNames[markIndex], assayStart = assayStart));
	} else {
		prompt = if (is.null(assayStart)) {""}
		else if (sum(markNames %in% assayStart) == 0) {paste('  No default assay start marks found.\n', sep = "")}
		else {paste('  Two or more default mark names found.\n')}
		prompt = paste(prompt, '  Mark names found:\n  "', paste(markNames, collapse = '" "'), '"\n', sep = "");
		cat(prompt);
		prompt = '  Which mark is the assay start? (enter "q" to skip assay start for this log or press ESC to abort)\n  > ';
		userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
		if (userInput == "q") return(list(time = 0, name = NA, assayStart = assayStart));
		userInput = .autocomplete(userInput, markNames);
		while (!(userInput %in% markNames)) {
			userInput = gsub('^["\']','', gsub('["\']$','', readline('  Please enter a valid mark name or "q", or press ESC to abort: ')));
			if (userInput == "q") return(list(time = 0, name = NA, assayStart = assayStart));
			userInput = .autocomplete(userInput, markNames);
		}
		
		if (!userInput %in% assayStart) {
			addMark = .getYesOrNo(paste('  Do you want to save "', userInput, '" as a default assay start mark? ', sep = ""))
			if (addMark) assayStart = c(assayStart, userInput);
		}
		
		markIndex = which(markNames == userInput);
		return(list(time = .timeToSeconds(marks[[markIndex]][2]), name = markNames[markIndex], assayStart = assayStart));
	}
}

# Prints out .findDupBehaviors in a more readable format, preceded by two spaces.
.printFindDupBehaviors = function(data) {
	behTable = .findDupBehaviors(data);
	for (i in 1:length(behTable)) {
		cat("  \"", names(behTable)[i], "\" (occurs in ", behTable[i], " logs)\n", sep = "");
	}
}

# Helper function for .getDataBatch()
# Provides a nice user interface for merging duplicate behaviors. Returns <data> with the
# dupped behaviors merged.
.promptToElimDups = function(data) {
	prompt = 'Please enter the name of a behavior that you would like to merge into another behavior, or "l" to print the current list of behaviors or "s" to save and quit.\n> '
	userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
	userInput = .autocomplete(userInput, c(names(.findDupBehaviors(data)), "s", "l"), caseSensitive = TRUE);
	while (userInput != "s") {
		if (userInput == "l") {
			.printFindDupBehaviors(data);
		} else if (userInput %in% names(.findDupBehaviors(data))) {
			prompt2 = paste('  What behavior do you want to merge "', userInput, '" with? (enter "q" to cancel)\n  > ', sep = "");
			input2 = gsub('^["\']','', gsub('["\']$','', readline(prompt2)));
			input2 = .autocomplete(input2, c(names(.findDupBehaviors(data)), "q"), caseSensitive = TRUE);
			while (!(input2 %in% c(names(.findDupBehaviors(data)), "q"))) {
				input2 = gsub('^["\']','', gsub('["\']$','', readline(prompt2)));
				input2 = .autocomplete(input2, c(names(.findDupBehaviors(data)), "q"), caseSensitive = TRUE);
			}
			confirmPrompt = paste('  Are you sure you want to merge behavior "', userInput, '" into "', input2,'"? ', sep = "")
			if (input2 != "q" && .getYesOrNo(confirmPrompt)) {
				data = .replaceBehAll(data, userInput, input2);
			}
		} else {
			cat("  That is not a valid behavior name.\n")
		}

		userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
		userInput = .autocomplete(userInput, c(names(.findDupBehaviors(data)), "s", "l"), caseSensitive = TRUE);
	}
	data = .filterDataList(data, renameSubjects = T);
	return(data);
}

# Source: ethograms_from_scorevideo.R
# Helper function for .getData
# Gets description table. Currently not called by anything.
.getDescriptionTableFromRawData = function (data0) {
	desc_row = which(grepl('^start.*description$', data0[,1])) + 2;
	desc_table = matrix(nrow=1, ncol=3);
	while (TRUE) {
		if (grepl('^-.*-$', data0[desc_row, ])) {
			break;
		} else {
			tmp = unlist(strsplit(data0[desc_row, ], '[0123]'));
			code = unlist(strsplit(gsub(' ','',tmp[1]),''));
			if (length(code) > 1) {
				code <- code [1];
			}
			description = gsub('  ','',tmp[length(tmp)]);
			subject = unlist(strsplit(data0[desc_row, ], '[^0123]'));
			desc_table = rbind(desc_table, c(code, description, subject[15]));
			desc_row = desc_row + 1;
		}
	}
	desc_table = desc_table[-1, ];
	dimnames(desc_table)[[2]] <- c("code", "description", "subj")
	return(desc_table);
}

# Helper function for .getData
# Parses full log to populate data frame
.parseFullLog = function (data0, desc_table, framesPerSecond) {
	start_full = which(data0=='FULL LOG') + 4;
	end_full = which(data0=='NOTES') - 2;
	if (start_full > end_full) return(NA);

	full = strsplit(data0[start_full:end_full, ], '    ');
	full = lapply(full, function(f) gsub('^ *','', gsub(' *$','', f)));
	full = lapply(full, function(f) f[f!=""]);
	
	# Get time, behavior, subject, type
	df = data.frame(matrix(nrow=1, ncol=6));
	for (entry in full) {
		isStart = if (length(entry)==5) {entry[5]} else {NA}
		newRow = c(as.numeric(entry[1]), as.character(entry[3]), as.character(entry[4]), isStart, NA, NA)
		df = rbind(df, newRow);
	}
	names(df) = c('time', 'behavior', 'subject', 'type', 'pair_time', 'duration');
	df$type <- as.character(df$type);
	df$time <- as.numeric(df$time) / framesPerSecond;
	df$pair_time <- as.numeric(df$pair_time);
	df$duration <- as.numeric(df$duration);
	df = df[-1, ];
	
	# Match up starts and stops; calculate durations
	for (i in 1:(dim(df)[1])) {
		entry = df[i,];
		if (!is.na(entry$type) && entry$type == 'stop') {
			this_time = entry$time;
			startIndex = 0;
			for (j in i:1) {
				if (df[j,]$behavior == entry$behavior && df[j,]$type=='start') {
					startIndex = j;
					break;
				}
			}
			df$pair_time[startIndex] <- entry$time;
			df$duration[startIndex] <- df$pair_time[startIndex] - df$time[startIndex]
			df$pair_time[i] <- df$time[startIndex]
			df$duration[i] <- df$duration[startIndex]
		}
	}
	return(df);
}

# Provides an interactive interface to stitch logs together that come from the same
# subject. Uses regexes to identify subject - the regex used becomes the name of the
# log. This allows conservation of folder structure!
# If logOutPref is provided, the stitched-together logs are written out as tables.
.stitchLogsTogether = function(data, logOutPref = NULL) {
	cat("Log names:\n");
	print(names(data));
	
	prompt = "  Please enter a subject name (regex), or \"l\" to print the current list of logs or \"q\" to quit logStitcher.\n  > ";
	regex = readline(prompt);
	
	while (regex != "q") {
		if (regex == "l") {
			cat("Log names:\n")
			print(names(data));
		} else if (regex != "") {
			matchingNames = which(grepl(regex, names(data)));
			if (length(matchingNames) > 0) {
				cat("  Logs that belong to subject ", regex, ":\n", sep = '"');
				print(names(data)[matchingNames]);
				if (.getYesOrNo("  Is this correct? ")) {
					if (!.logsOverlap(data[matchingNames]) || .getYesOrNo("Proceed anyway? ")) {
						newSubDat = data[[matchingNames[1]]];
						if (length(matchingNames) > 1) {
							for (i in 2:length(matchingNames)) newSubDat = rbind(newSubDat, data[[matchingNames[i]]]);
						}
						
						newSubDat = .sortByTime(newSubDat);
						data = data[-matchingNames];
						data[[length(data) + 1]] <- newSubDat;
						names(data)[length(data)] <- regex;
					}
				}
			} else {
				cat("  No logs belong to subject ", regex, ".\n", sep = '"');
			}
		}
		regex = readline(prompt);
	}
	
	if (!is.null(logOutPref)) {
		for (i in 1:length(data)) {
			write.table(data[[i]], file = paste(logOutPref, names(data)[i], sep = ""), sep = "\t");
		}
	}
	return(data);
}

# Helper function for .stitchLogsTogether
# Returns FALSE if no logs in <logs> have overlapping times, and TRUE if they do.
# Additionally throws a warning for each set of overlapping logs found.
.logsOverlap = function(logs) {
	if (length(logs) < 2) return(FALSE);
	logMins = unlist(lapply(logs, function(d){min(d$time)}));
	logMaxes = unlist(lapply(logs, function(d){max(d$time)}));
	overlapFound = FALSE;
	for (i in 1:(length(logs) - 1)) {
		for (j in (i+1):length(logs)) {
			if (logMaxes[i] > logMins[j] && logMins[i] < logMaxes[j]) {
				warning('Logs "', names(logs)[i], '" and "', names(logs)[j], '" overlap.', immediate. = TRUE)
				overlapFound = TRUE;
			}
		}
	}
	return(overlapFound);
}


.pairLogs = function(data) {
	groupwiseData = .sepGroups(data);
	groupNames = groupwiseData$groupNames;
	groupData = groupwiseData$groupData;
	if (length(groupNames) != 2) stop("This code is written to deal with two groups. Please see Katrina.");
#    $groupNames, the names of the groups;
#    $groupData, a list of lists of data frames. groupData[[1]] is a list of the data frames of all the
#      subjects in group 1, groupData[[2]] has the data for the subjects in group 2, etc.
#    $behnames, the behaviors that occur in any log.

	key = .makePairKey(names(groupData[[1]]), names(groupData[[2]]), groupNames);
	return(data[as.vector(key)]);
}

.makePairKey = function(group1Names, group2Names, groupNames) {
	key = matrix(ncol = 2, nrow = length(group1Names));
	dimnames(key)[[2]] <- groupNames;
	key[,1] <- group1Names;
	group1Names = gsub(paste("^", groupNames[1], "/", sep = ""), "", group1Names);
	# print(key);
	for (i in 1:length(group1Names)) {
		subLen = 5;
		uniqueMatches = NULL;
		while (is.null(uniqueMatches)) {
			for (j in 1:(nchar(group1Names[i]) - subLen)) {
				searchstring = substr(group1Names[i], j, j + subLen);
				matches = grepl(paste("^", groupNames[2], "/.*", searchstring, ".*$", sep = ""), group2Names);
				if(sum(matches) == 1) {
					uniqueMatches = c(uniqueMatches, group2Names[matches]);
					uniqueMatches = names(table(uniqueMatches));
				}
			}
			if (!is.null(uniqueMatches) && length(uniqueMatches) > 1) uniqueMatches = NULL;
			subLen = subLen + 1;
		}
		
		key[i,2] <- uniqueMatches;
	}
	return(key);
}


#####################################################################################################
## FILTERING, SORTING, AND EDITING DATA                                                            ##
#####################################################################################################

######################################### SORTING ###################################################

# Returns <data> sorted in order of increasing <attribute>. <attribute> should be a vector
# of numbers, strings, or logical values, the same length as <data> with value attribute[i]
# corresponding to log data[[i]].
# Arguments can be passed through the ... to order(), including (notably) decreasing=T which
# sorts in order of decreasing <attribute> and na.last = NA which removes NAs from the data.
# Or na.last = T to put NAs last, or na.last = F to put them first.
# TODO error checking
.sortByAttribute = function(data, attribute, ...) {
	if (is.character(attribute)) attribute = as.factor(attribute);
	# print(attribute[order(attribute, ...)]);
	return(data[order(attribute, ...)]);
}

# Extracts an attribute vector to use in sorting from a matrix or data frame of
# supplemental data (<suppData>). <data> should be the list of score logs. You must
# also provide <attributeCol> which is either the column index or the column name of
# suppData that contains the data of interest (for example GSI), and either <nameCol>
# (the name or index of a column in suppData with the names of the score logs) or
# <indexCol> (the name or index of a column in suppData with the index in <data> of
# each row's associated score log).
# TODO NOTE: In helpfile, mention that suppData needs to be read with read.csv.
.getAttributeFromSupplementalData = function(suppData, data, attributeCol, nameCol = NA, indexCol = NA) {
	if (!.checkColumnOK(attributeCol, suppData)) stop("Invalid attributeCol.");
	if (is.na(indexCol)) {
		if (is.na(nameCol)) stop("You must provide either nameCol or indexCol.");
		if (!.checkColumnOK(nameCol, suppData)) stop("Invalid nameCol.");
		suppData = as.data.frame(suppData);
		for (i in 1:length(suppData[,nameCol])) {
			subj = suppData[i, nameCol];
			logs = grepl(subj, names(data));
			if (sum(logs) == 0) suppData$index.tmp[i] <- NA
			else if (sum(logs) == 1) suppData$index.tmp[i] <- which(logs)
			else warning(paste("More than one log matches assay name", subj), immediate. = TRUE)
			indexCol = which(names(suppData) == "index.tmp");
		}
	} else {
		if (!.checkColumnOK(indexCol, suppData)) stop("Invalid indexCol.")
		suppData[,indexCol] <- as.numeric(suppData[,indexCol]);
	}
	
	attribute = NULL;
	for (i in 1:length(data)) {
		ai = suppData[,attributeCol][which(suppData[,indexCol] == i)];
		if (length(ai) != 1) {
			warning("No row found in suppData for log ", names(data)[i], immediate.=TRUE);
			ai = NA;
		}
		attribute = c(attribute, ai);
	}
	return(attribute);
}

# The next four functions return behavior count, latency, total duration, and average duration
# as an attribute-vector to use in .sortByAttribute().

.getCountAttribute = function(behavior, data, warn = T) {
	counts = rep(NA, length(data));
	for (i in 1:length(data)) {
		counts[i] = sum(data[[i]]$behavior == behavior);
	}
	numNot0s = sum(counts != 0);
	if (numNot0s == 0) stop("No occurances of \"", behavior, "\" found.")
	else if (numNot0s < length(data) / 2 && warn) warning("Behavior \"", behavior, "\" only found in ", numNot0s, " out of ", length(data), " logs.", immediate.=TRUE);
	return(counts);
}

.getLatencyAttribute = function(behavior, data, n = 1) {
	latencies = rep(NA, length(data));
	for (i in 1:length(data)) {
		behOccurances = data[[i]]$time[data[[i]]$behavior == behavior];
		if (length(behOccurances) >= n) latencies[i] <- behOccurances[n];
	}
	numNotNAs = sum(!is.na(latencies));
	if (numNotNAs == 0) stop("No occurances of \"", behavior, "\" found.")
	else if (numNotNAs < length(data) / 2) warning("Behavior \"", behavior, "\" only found in ", numNotNAs, " out of ", length(data), " logs.", immediate.=TRUE);
	return(latencies);
}

.getTotalDurAttribute = function(behavior, data) {
	durations = rep(0, length(data));
	for (i in 1:length(data)) {
		behOccurances = data[[i]][data[[i]]$behavior == behavior,];
		if (length(behOccurances$duration) > 0) {
			if (sum(is.na(behOccurances$type)) > 0) stop("Behavior \"", behavior, "\" is not durational.");
			durations[i] = sum(behOccurances$duration[which(behOccurances$type == "start")]);
		}
	}
	numNot0s = sum(durations != 0);
	if (numNot0s == 0) stop("No occurances of \"", behavior, "\" found.")
	else if (numNot0s < length(data) / 2) warning("Behavior \"", behavior, "\" only found in ", numNot0s, " out of ", length(data), " logs.", immediate.=TRUE);
	return(durations);
}

.getAverageDurAttribute = function(behavior, data) {
	totalDurations = .getTotalDurAttribute(behavior, data);
	counts = .getCountAttribute(behavior, data, warn = F);
	zeroCounts = which(counts == 0);
	if (sum(totalDurations[zeroCounts]) != 0) stop("NONZERO DURATION WITH ZERO COUNT. This is a bug.");
	counts[zeroCounts] <- 1; # to avoid div-by-0
	return(totalDurations / counts);
}


######################################## FILTERING ##################################################

# OPTIONS:
# startTime, endTime - only get behavior from times [startTime, endTime]
# subjects - only consider these subjects' behavior
# startOnly - if TRUE, only consider starts of behaviors (ignore ends). If a character vector,
#   ignore ends only for the behaviors named in the character vector.
# boutInterval - interval to separate bouts, in seconds. Default is null (no bout separation)
# stateBehaviors - when boutInterval is given, a vector of behaviors to treat as states. These
#   behaviors can start in one bout and end in another bout.
# minNumBehaviors - remove behaviors that do not occur at least this many times.
# toExclude - remove all behaviors in this character vector.
# splitPot - option specifically for Scott's PGF2a data to separate "inside POT" by subject into
#   into "male IN POT" and "female IN POT"
# renameStartStop - for durational behaviors, append "start" to the behavioral description at
#   the start of the behavior, and "stop" to the behavioral description at the end of the behavior.
# zeroBeh - Changes time=0 in all logs from being the start of the assay to being the first
#   of <zeroBeh>. If <zeroBeh> never occurs, the whole log is in negative time. Runs AFTER
#   startTime/endTime!
.filterData = function(data, startTime = NA, endTime = NA, subjects = NULL, startOnly = NULL,
					   boutInterval = NULL, stateBehaviors = NULL, minNumBehaviors = NULL, toExclude = NULL,
					   renameStartStop = FALSE, zeroBeh = NULL, zeroBehN = 1) {
	if (!is.na(startTime)) {data <- data[data$time >= startTime,];}
	if (!is.na(endTime)) {data <- data[data$time <= endTime,];}
	if (!is.null(subjects)) {
		data <- data[data$subject %in% subjects,]
		if(length(data$behavior) == 0) {stop(paste("Error: Incorrect Subject Name:", subjects))}
	}
	if (!is.null(startOnly)) {
		if (is.logical(startOnly) && startOnly) {
			data <- data[is.na(data$type) | data$type != "stop",];
			data$type <- NA;
		} else if (!is.logical(startOnly)) {
			for (beh in startOnly) {
				data$type[data$behavior == beh & !is.na(data$type) & data$type == "start"] <- NA;
				data <- data[!(data$behavior == beh & !is.na(data$type) & data$type == "stop"),]
			}
		}
	}
	if (!is.null(boutInterval)) {
		data = .separateBouts(data, boutInterval);
	}
	if (!is.null(minNumBehaviors)) {
		freqtable = table(data$behavior);
		behaviorsToTrash = names(freqtable)[freqtable < minNumBehaviors]
		data <- data[!(data$behavior %in% behaviorsToTrash),];
		if(length(data$behavior) == 0) {stop(paste("Error: minNumBehaviors is set too high. No behavior occurs", minNumBehaviors, "times."))}
	}
	if (!is.null(toExclude)) {
		data <- data[!(data$behavior %in% toExclude),];
	}
	if (renameStartStop) {data = .renameStartStop(data);}
	if (!is.null(zeroBeh)) {data = .setZeroToNthOfBehavior(data, zeroBeh, zeroBehN)}
	return(data);
}

# lapply's .filterData with the selected options to a list of score logs.
# The option renameSubjects searches all the logs for any two behaviors that have identical names
#   but different subjects. If any such behaviors are found, they are replaced by "subject old_behavior_name"
#   to disambiguate. 
# renameSubjects is here because we want to rename the same behaviors in the same way across all logs.
.filterDataList = function(data, renameSubjects = F, ...) {
	if (renameSubjects) {
		behaviors <- names(.findDupBehaviors(data));
		problemBehs = character()
		for (beh in behaviors) {
			if (length(table(unlist(lapply(data, function(d) {d[d$behavior == beh,]$subject})))) != 1) problemBehs <- c(problemBehs, beh);
		}
		if (length(problemBehs > 0)) {
			cat("Renaming behaviors: \"");
			cat(problemBehs, sep = '" "');
			cat('"\n')
		}
		data = lapply(data, function(d) {.sepSubject(d, problemBehs)});
	}
	return(lapply(data, function(d) {.filterData(d, ...)}));
}

# Removes empty score logs from <data>. BE CAREFUL WITH THIS - it might mess up all kinds
# of things by changing the number and order of logs.
.removeNALogs = function(data) {
	emptyLogs = which(unlist(lapply(data, function(d){!is.data.frame(d)})))
	return(data[-emptyLogs]);
}

# Might someday become a handy interface for .filterDataList(). Or nah. TODO complete or trash
# .interactiveFilter = function(data) {
	# allowableCommands = c("help", "time binning", "ignore behaviors by subject(s)", "make non-durational", "separate bouts", "rename behavior",
						  # "ignore infrequent behaviors", "ignore behavior(s)", "center logs around nth occurance of a behavior");

	# cat("Welcome to interactive data filter. Please type a command to begin, or \"help\" to get a list of commands.\n> ");
# }


# Renames the behaviors in <toSeparate>, which occur in two or more subjects, so that the behavior
# description indicates which subject it is.
.sepSubject = function(data, toSeparate) {
	toSep = data$behavior %in% toSeparate;
	data$behavior[toSep] <- paste(as.character(data$subject[toSep]), data$behavior[toSep]);
	return(data);
}

# Renames behaviors to differentiate between starts and stops. If the behavior has already been renamed,
# this function does nothing.
.renameStartStop = function(data) {
	toReplace = !is.na(data$type) & !grepl(" st[oa][rp]t?$", data$behavior);
	data$behavior[toReplace] <- paste(data$behavior[toReplace], data$type[toReplace]);
	return(data);
}

# Finds the nth occurance of <behavior> in the log <data> and sets data's <time> column to the offset
# relative to that behavior.
.setZeroToNthOfBehavior = function(data, behavior, n = 1) {
	data = rbind(data, list(time = 0, behavior = "assay start", subject = NA, type = NA, pair_time = NA, duration = NA)); # TODO make this line better
	targetBehIndices = which(data$behavior == behavior);
	newStartTime = 0;
	if (length(targetBehIndices) >= n) {
		newStartTime = data$time[targetBehIndices[n]];
	} else {
		newStartTime = max(data$time) + 1;
	}
	data$time <- data$time - newStartTime; # TODO change attribute
	return(data);
}


# Whenever there are <intervalToSeparate> seconds between adjacent behaviors, inserts a "STOP" after the
# behavior before the pause and a "START" before the behavior after the pause. Also inserts a "START" at the
# beginning of the log and a "STOP" at the end.
# TODO use type, pair_time, duration to pair "START"s with "STOP"s.
.separateBouts = function (data, intervalToSeparate, stateBehaviors = NULL) {
	# 	names(df) = c('time', 'behavior', 'subject', 'type', 'pair_time', 'duration');
	data$type[!is.na(data$type) & data$type == "start" & data$behavior %in% stateBehaviors] <- "sTaRt";
	data$type[!is.na(data$type) & data$type == "stop" & data$behavior %in% stateBehaviors] <- "sToP";
	
	newData = data.frame(time = data$time[1], behavior = "START", subject = NA, type = NA, pair_time = NA, duration = NA); 
	newData = rbind(newData, data[1,]);
	# print(data[1:20,]);
	for (i in 2:length(data[,1])) {
	#	print(paste(sum(data$type[!is.na(data$type)][1:(i-1)] == "start"), sum(data$type[!is.na(data$type)][1:(i-1)] == "stop"), sum(data$type[!is.na(data$type)][1:(i-1)] == "start") == sum(data$type[!is.na(data$type)][1:(i-1)] == "stop")));
		if (as.numeric(data$time[i]) - as.numeric(data$time[i-1]) >= intervalToSeparate &&
		    sum(data$type[1:(i-1)][!is.na(data$type)[1:(i-1)]] == "start") == sum(data$type[1:(i-1)][!is.na(data$type)[1:(i-1)]] == "stop")) {
			stopRow = data.frame(time = data$time[i-1], behavior = "STOP", subject = NA, type = NA, pair_time = NA, duration = NA);
			newData = rbind(newData, stopRow);
			startRow = data.frame(time = data$time[i], behavior = "START", subject = NA, type = NA, pair_time = NA, duration = NA);
			newData = rbind(newData, startRow);
		}
		newData = rbind(newData, data[i,]);
	}
	stopRow = data.frame(time = data$time[length(data$time)], behavior = "STOP", subject = NA, type = NA, pair_time = NA, duration = NA);
	newData = rbind(newData, stopRow);
	dimnames(newData)[[1]] <- 1:length(dimnames(newData)[[1]])
	
	newData$type[!is.na(newData$type) & newData$type == "sTaRt"] <- "start";
	newData$type[!is.na(newData$type) & newData$type == "sToP"] <- "stop";
	return(newData);
}

# Makes <leader> into a durational behavior with stops at each occurance of <follower> in the log <data>.
# If <rename>, follower occurances are renamed to <leader>. (recommended)
.makeDurationalBehavior = function (data, leader, follower, rename = T) {
	if (leader == follower) stop("Leader cannot be the same as follower.");
	leaderIndices = which(data$behavior == leader);
	followerIndices = which(data$behavior == follower);
	
	data$type[leaderIndices] <- "start";
	data$type[followerIndices] <- "stop";
	if (rename) data$behavior[followerIndices] <- leader;
	
	if (length(leaderIndices) == 0) {
		if (length(followerIndices) > 0) {
			warning("Follower occurs but leader doesn't.", immediate. = T)
		}
		return(data);
	}
	
	if (length(followerIndices) >= 1 && followerIndices[1] < leaderIndices[1]) {
		if (sum(followerIndices < leaderIndices[1]) == 1) {
			warning("Follower at the beginning of log without a leader.", immediate. = T);
		} else {
			stop("Two occurances of follower before first occurance of leader.");
		}
	}
	for (i in 1:length(leaderIndices)) {
		thisLeader = leaderIndices[i];
		nextLeader = if (i == length(leaderIndices)) length(data$behavior) + 1 else leaderIndices[i+1];
		followers = followerIndices[followerIndices > thisLeader & followerIndices < nextLeader];
		if (length(followers) < 1) {
			if (i == length(leaderIndices)) {
				warning("Leader at the end of log without a follower.", immediate. = T);
			} else {
				stop("Two occurances of leader in a row (indices ", thisLeader, " and ", nextLeader, ")");
			}
		} else if (length(followers) > 1) {
			stop("Two occurances of follower in a row (indices ", followers[1], " and ", followers[2], ")");
		} else {
			data$pair_time[thisLeader] <- data$time[followers];
			data$pair_time[followers] <- data$time[thisLeader];
			data$duration[c(thisLeader, followers)] <- data$time[followers] - data$time[thisLeader];
		}
	}

	return(data);
}




######################################### EDITING ###################################################

# Replaces behavior description <toReplace> with description <replacement> in data
# frame data. An example use of this function would be to replace "Male in pot" with
# "Male In Pot"
.replaceBeh = function(data, toReplace, replacement) {
	data$behavior[data$behavior == toReplace] <- replacement;
	return(data);
}

# Calls .replaceBeh on every data frame of a data list.
.replaceBehAll = function(data, toReplace, replacement) {
	for (i in 1:length(data)) {
		data[[i]] <- .replaceBeh(data[[i]], toReplace, replacement);
	}
	return(data);
}

# Returns a table with names of all the behaviors in the logs in <data>, and counts of
# how many logs each behavior appears in.
.findDupBehaviors = function(data) {
	return(table(unlist(lapply(data, function(f) {names(table(f$behavior))}))));
}


# This is an example of how to use .replaceBehAll to clean up data such that all logs
# have the same descriptions for the same behavior. This accomplishes the task for the
# logs provided by Mariana.
# To find behaviors in your data that may be named differently in different logs, call
# .findDupBehaviors. You can visually inspect this output to look for duplicates.
.cleanUpMariana = function(mariana_dat) {
	mariana_clean <- .replaceBehAll(mariana_dat, "Male darts", "Male Darts");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male Darting", "Male Darts");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male in Pot", "Male In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male in pot", "Male In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male to Pot after Female", "Male In Pot")
	mariana_clean <- .replaceBehAll(mariana_clean, "Female flees", "Female Flees");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female follows", "Female Follows");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female in pot", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female in Pot", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female to Pot w/o Male", "Female In Pot");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male bites/rams", "Male Bites/Rams");
	mariana_clean <- .replaceBehAll(mariana_clean, "Male chases", "Male Chases");
	mariana_clean <- .replaceBehAll(mariana_clean, "Female Acceptance Still, Faces Male", "Female Follows");
	return(mariana_clean);
}


# Cleans up Scott's pgf2a data.
.cleanUpPGF2a = function(my_data) {
	my_data <- .filterDataList(my_data, toExclude=c("approach", "APPROACH", "BITE", "flee", "QUIVER"))
	my_data <- .filterDataList(my_data, startOnly=c("chase", "female FOLLOW", "FLEE", "FOLLOW", "lead", "male CHASE", "male LEAD", "male QUIVER", "quiver"))
	.cleanFxn = function(dat) {
		clean <- .replaceBehAll(dat, "female FOLLOW", "female follow");
		clean <- .replaceBehAll(clean, "female IN POT", "female in pot");
		clean <- .replaceBehAll(clean, "FLEE", "female flee");
		clean <- .replaceBehAll(clean, "FOLLOW", "female follow");
		clean <- .replaceBehAll(clean, "inside pot", "male in pot");
		clean <- .replaceBehAll(clean, "inside POT", "female in pot");
		clean <- .replaceBehAll(clean, "male BITES", "bite");
		clean <- .replaceBehAll(clean, "male CHASE", "chase");
		clean <- .replaceBehAll(clean, "male IN POT", "male in pot");
		clean <- .replaceBehAll(clean, "male LEAD", "lead");
		clean <- .replaceBehAll(clean, "male QUIVER", "quiver");
		return(clean);
	}
	my_data <- .cleanFxn(my_data);
	return(my_data);
}


#####################################################################################################
## GETTING STATS                                                                                   ##
#####################################################################################################

# Takes a list of data frames and sorts them into groups based on the filenames. The name of each
# group will be the name of the folder the score logs came from in .getDataBatch(). Returns a list
# containing:
#    $groupNames, the names of the groups;
#    $groupData, a list of lists of data frames. groupData[[1]] is a list of the data frames of all the
#      subjects in group 1, groupData[[2]] has the data for the subjects in group 2, etc.
#    $behnames, the behaviors that occur in any log.
# In this implementation, score logs that did not come from a folder (either from grouped datasets
#    or ungrouped datasets) are placed in "Default Group".
.sepGroups = function(data) {
	groupNames = names(table(gsub("((.+)/)?.+", "\\2", names(data))));
	ungrouped = "" %in% groupNames
	groupNames[groupNames == ""] <- "Default Group";
	behnames = names(table(unlist(lapply(data, function(f) {names(table(f$behavior))}))));
	groupData = list();
	for (i in 1:length(groupNames)) {
		groupData[[i]] = data[grepl(paste("^", groupNames[i], sep = ""), names(data))];
		names(groupData)[i] <- groupNames[i];
	}
	if (ungrouped) {
		defaultGroupIndex = groupNames[groupNames == "Default Group"];
		groupData[[defaultGroupIndex]] = data[grepl("^[^/]*$", names(data))];
	}
	return(list(groupNames = groupNames, groupData = groupData, behnames = behnames));
}

# Helper function for .calcBasicStats(). Loops through every data frame in the list <data> and gets the
# count and latency of each behavior in <behnames>, and the total duration of each behavior in <durBehNames>.
# If a behavior does not occur in a score log, the latency will be NA.
# Returns a list containing a matrix of each type of data, as well as a matrix containing all the types of
# data, sorted by behavior (ready to be written to a .csv)
.extractBasicStats = function(data, behnames, durBehNames) {
	startOnlyDat = .filterDataList(data, startOnly = T);
	tables = lapply(startOnlyDat, function(f) {table(f$behavior)});
	countsMat = matrix(nrow = length(names(tables)), ncol = length(behnames), dimnames = list(names(tables), paste(behnames, ": Count")));
	latMat = matrix(nrow = length(names(tables)), ncol = length(behnames), dimnames = list(names(tables), paste(behnames, ": Latency")));
	for (i in 1:length(names(tables))) {
		for (j in 1:length(behnames)) {
			if (behnames[j] %in% names(tables[[i]])) {
				countsMat[i,j] <- tables[[i]][[behnames[j]]];
				latMat[i,j] <- startOnlyDat[[i]][startOnlyDat[[i]]$behavior == behnames[j],][1,]$time;
			}
			else {
				countsMat[i,j] <- 0;
				#latMat[i,j] is already initialized to NA
			}
		}
	}
	
	durTMat = durAMat = NULL;
	if (!is.null(durBehNames)) {
		durTMat = matrix(nrow = length(names(tables)), ncol = length(durBehNames), dimnames = list(names(tables), paste(durBehNames, ": Duration Total")));
		durAMat = matrix(nrow = length(names(tables)), ncol = length(durBehNames), dimnames = list(names(tables), paste(durBehNames, ": Duration Average")));
		for (i in 1:length(names(tables))) {
			for (j in 1:length(durBehNames)) {
				behData = data[[i]][data[[i]]$behavior == durBehNames[j] & data[[i]]$type == "start",];
				durTMat[i,j] = sum(as.numeric(behData$duration));
				durAMat[i,j] = if(length(behData$behavior) == 0) NA else durTMat[i,j] / length(behData$behavior);
			}
		}
		durTMat = t(durTMat);
		durAMat = t(durAMat);
	}
	comboMat = rbind(t(countsMat), t(latMat), durTMat, durAMat);
	comboMat = comboMat[order(dimnames(comboMat)[[1]]),];
	
	return(list(counts = t(countsMat), latencies = t(latMat), durations = durTMat, total = comboMat));
}


# Calculates the count, latency, and duration for each behavior that occurs in any score log in <data>. Saves the
#   results to .csv files with prefix <outfilePrefix>. Also calculates the average and standard deviation for each
#   group for each measure of each behavior.
# This function also calculates p values comparing the first two groups (normally, there are only two groups - experimental
#   and control - so this behavior is not terrible) with three different tests - wilcox.test(), t.test(), and bootstrap2independent().
#   Latencies are only compared for behaviors that occur in at least <minNumLogs> score logs in each group. Graphs output by the
#   bootstrap function are also saved.
.calcBasicStats = function(data, outfilePrefix, tests = list(t.test = t.test, wilcox = wilcox.test, bootstrap = list(func = .bootstrapWrapper)), minNumLogs = 3) {
	groupwiseLogs = .sepGroups(data);

	durBehDat = lapply(data, function(d) {d[!is.na(d$type) & d$type == "start",]});
	durBehNames = names(table(unlist(lapply(durBehDat, function(f) {names(table(f$behavior))}))));
	dataByGroup = list();
	for (group in groupwiseLogs$groupNames) {
		dataByGroup[[group]] <- .extractBasicStats(groupwiseLogs$groupData[[group]], groupwiseLogs$behnames, durBehNames);
		write.csv(dataByGroup[[group]]$total, file = paste(outfilePrefix, group, "data.csv", sep = "_"));
	}
	return(.runStats(dataByGroup = lapply(dataByGroup, function(d){d$total}), outfilePrefix = paste(outfilePrefix, "basicstats", sep = "_"),
			tests = tests, minNumLogsForComparison = minNumLogs));
}

# A wrapper function for bootstrap2independent that makes it play well with .runStats
# argList must contain x, y, row, outfilePrefix, groupNames, and can optionally contain <trials> (the number of trials).
# The function must be passed in in a list.
.bootstrapWrapper = function(argList) {
	# print(argList);
	if(!("trials" %in% names(argList))) argList$trials <- 10000;
	bs = bootstrap2independent(argList$x, argList$y, dataDescriptor = argList$row,
	       						 outfile = paste(argList$outfilePrefix, gsub("[ :/]", "", argList$row), "bootstrap.jpg", sep = "_"),
	       						 groupNames = argList$groupNames, trials = argList$trials, printResults = FALSE, verbose = FALSE);
	# print(list(p = bs$p.value, dat = bs$data));
	return(list(p.value = bs$p));
}

.bootstrapPairedWrapper = function(argList) {
	if(!("trials" %in% names(argList))) argList$trials <- 10000;
	sumVec = argList$x - argList$y; # TODO BUG this should be addition (+). need to ask austin abt case (1,2,0) (1,2,0) first tho. pot exit --> flee from female
	if (sum(sumVec[!is.na(sumVec)]) == 0) return(list(p.value = NA));
	bs = bootstrap2paired(argList$x, argList$y, dataDescriptor = argList$row,
	       						 outfile = paste(argList$outfilePrefix, gsub("[ :/]", "", argList$row), "bootstrap.jpg", sep = "_"),
	       						 conditionNames = argList$groupNames, trials = argList$trials, printResults = FALSE, verbose = FALSE);
	# print(list(p = bs$p.value, dat = bs$data));
	return(list(p.value = bs$p));
}

# need to library(survival) TODO
# needs assayLength in arglist
.coxphWrapper = function(argList) {
	x = argList$x;
	y = argList$y;
	data = c(x,y);
	in.group1 = 1:length(data) <= length(x);
	censored = is.na(data);
	data[censored] <- argList$assayLength;
	survObj = Surv(data, !censored);
	testout <- coxph(survObj ~ in.group1);
	return(list(p.value = as.data.frame(coef(summary(testout)))$Pr));
}

# Calculates the average and standard deviation of <dataByGroup> (and optionally statistical tests in <tests>) and
#   outputs the results to a .csv starting with <outfilePrefix>. If <skipNA>, the average and standard deviation
#   are calculated ignoring NAs. ie, average = (sum of not-NA data) / (number of subjects with not-NA data)
# <dataByGroup> should be a list of matrices. Each matrix represents an experimental group, with a column for each
#   subject and a row for each measurement (ie "Lead count" or "Quiver -> Lead Transitional Probability"). The list
#   has length two for a two-group experiment.
# <tests> is a list of functions. The names of the list are used to label columns. e.g. list(ttest = t.text, wilcox = wilcox.test)
# The test functions must take their data as x and y and return a list with the p value stored under name "p.value".
# If a function needs other arguments than just the data:
#   Pass it in as a list(func = myFunctionWrapper, arg1Name = arg1, arg2Name = arg2, etc). The function MUST be the first argument.
#   The function will be called as myFunctionWrapper(list(arg1Name = arg1, arg2Name = arg2, etc, x = data, y = data, etc))
#   You also get the row name, group names, and outfile prefix passed along as part of that list. You can make a list of length
#     1 to trigger this behavior if that's all you needed anyway.
#   Inside myFunctionWrapper you might do something like:
#   myFunctionWrapper = function(arglist) {
#		return(myFunction(group1 = arglist$x, group2 = arglist$y, arg1Name = arglist$arg1Name, outpref = arglist$outfilePrefix))
#   }
.runStats = function(dataByGroup, outfilePrefix, tests, minNumLogsForComparison = 3, skipNA = F) {
	average = if (skipNA) {lapply(dataByGroup, apply, 1, function(row){mean(row[!is.na(row)])})} else lapply(dataByGroup, apply, 1, mean);
	stddev = if (skipNA) {lapply(dataByGroup, apply, 1, function(row){sd(row[!is.na(row)])})} else lapply(dataByGroup, apply, 1, sd);
	if (skipNA) {
		for (i in 1:length(dataByGroup)) {
			allNARows = apply(dataByGroup[[i]], 1, function(row){sum(!is.na(row)) == 0});
			average[[i]][allNARows] <- NA;
			stddev[[i]][allNARows] <- NA;
		}
	}
	rownames = dimnames(dataByGroup[[1]])[[1]];
	df = data.frame(average = average, stddev = stddev);
	offset = dim(df)[2];
	
	if (length(dataByGroup) >= 2 && length(tests) > 0) {
		if (length(dataByGroup) > 2) warning(paste("It looks like you have more than two experimental groups. ",
												   "Contact Katrina if you want code to deal with that.\nFor now, I'll run tests on the first two groups (\"",
												   names(dataByGroup)[1], "\" and \"", names(dataByGroup)[2], "\")", sep = ""), immediate. = TRUE)
		for (i in 1:length(tests)) {
			df = data.frame(df, numeric(dim(df)[1]));
			names(df)[i + offset] <- names(tests)[i];
		}
		for (i in 1:length(rownames)) {
			row = rownames[i]
			cat('Analyzing "', row, '"...\n', sep = "");
			group1dat = dataByGroup[[1]][row,];
			group2dat = dataByGroup[[2]][row,];
			
			if (sum(group1dat[!is.na(group1dat)]) + sum(group2dat[!is.na(group2dat)]) == 0) {
				df[i, (offset + 1):(offset + length(tests))] <- NA;
				next;
			}
			
			enoughObservations = (sum(!is.na(group1dat)) >= minNumLogsForComparison && sum(!is.na(group2dat)) >= minNumLogsForComparison);
			if (enoughObservations) {
				for (j in 1:length(tests)) {
					if (!is.list(tests[[j]])) {
						fxn = tests[[j]];
						df[i, j + offset] <- fxn(x=group1dat, y=group2dat)$p.value
					} else {
						functionList = tests[[j]];
						df[i, j + offset] <- functionList[[1]](c(functionList[-1], list(x = group1dat, y = group2dat, row = row,
															   outfilePrefix = outfilePrefix, groupNames = names(dataByGroup)[1:2])))$p.value
					}
				}
			} else {
				warning(paste('Skipping "', row, '" (not enough observations)', sep = ""), immediate.=T);
				df[i, (offset + 1):(offset + length(tests))] <- NA;
			}
		}
	}
	write.csv(df, file = paste(outfilePrefix, "stats.csv", sep = "_"));
	return(df);
}




#####################################################################################################
## PROBABILITY MATRICES                                                                            ##
#####################################################################################################

# Source: ethograms_from_scorevideo.R
# Reads in a vector giving a sequence of behaviors and returns a matrix giving the transitional
#  probability for each pair of behaviors.
# Usage: .getProbabilityMatrix(.renameStartStop(dataFrame)$behavior)      to include behavior starts and ends
.getProbabilityMatrix = function (data, removeZeroCol=F, ...) {
	beh = names(table(data[-length(data)]));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeTransitionProbability(data=data, leader=leader, follower=follower, ...);
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$probability;
		}
	}
	if (removeZeroCol) {
		colSums = apply(probMat, 2, sum);
		if (any(colSums==0)) {
			probMat = probMat[-which(colSums==0), -which(colSums==0)];
		}
	}
	return(probMat);
}

# Source: ethograms_from_scorevideo.R
# Helper function for .getProbabilityMatrix
# Computes the transition probability between leader and follower
.computeTransitionProbability = function (data, leader, follower, byTotal = FALSE) {
	count = 0;
	termination = 0;
	for (i in 1:length(data)) {
		if (data[i] == leader) {
			if (i == length(data)) {
				termination = 1;
			} else if (data[i+1] == follower) {
				count = count + 1;
			}
		}
	}
	total_leader = if(termination) {sum(data == leader) - 1} else {sum(data == leader)};
	prob = if (byTotal) {count / (length(data) - 1)} else if(total_leader) {count / total_leader} else {0} ;
	return(list(probability=prob, termination=termination, count_transitions=count, count_leader=total_leader));
}

# Calls .combineProbabilityMatrices on each experimental group in <data>, and returns a list of the group-level
# probability matrices for each group. If <byTotal>, the probability is
#	(# of occurances of transition from leader to follower) / (num transitions) ;
# otherwise, (default) the probability is
#	(# of occurances of transition from leader to follower) / (num occurances of leader).
# Notice that the combined probability is (probability for 1st subj + probability for 2nd subj + ...) / (number of subjects), NOT
# (total # of transitions from leader to follower in all logs) / (total (# of occurrances of leader or # of transitions) in all logs).
.groupLevelProbMats = function(data, byTotal = FALSE) {
	data = .filterDataList(data, renameStartStop=TRUE);
	groupwiseLogs = .sepGroups(data);
	# 	return(list(groupNames = groupNames, groupData = groupData, behnames = behnames));
	
	probMatsByGroup = list();
	for (group in groupwiseLogs$groupNames) {
		groupPMsAndCounts = list(probMats = lapply(groupwiseLogs$groupData[[group]], function(d) {.getProbabilityMatrix(d$behavior, byTotal=byTotal)}),
								 counts = lapply(groupwiseLogs$groupData[[group]], function(d) {table(d$behavior)}));
		probMatsByGroup[[group]] <- .combineProbabilityMatrices(groupPMsAndCounts, groupwiseLogs$behnames, byTotal);
	}
	return(probMatsByGroup);
}

# Helper function for .groupLevelProbMats()
# Takes <data>, a list of (1) a list of the probability matrices for each subject and (2) a list of the counts of
# each behavior for each subject, and combines them into a unified probability matrix.
.combineProbabilityMatrices = function(data, behnames = NULL, byTotal) {
	if (is.null(behnames)) {behnames = names(table(unlist(lapply(probmas, function(d){dimnames(d)[[1]]}))));}
	probmas = data$probMats;
	counts = data$counts;
	probMat = matrix(data = numeric(length(behnames) * length(behnames)),
					 nrow = length(behnames), ncol = length(behnames), dimnames = list(behnames, behnames));
	countVec = numeric(length(behnames));
	names(countVec) <- behnames;
	nfish <- countVec;
	for (row in rownames(probMat)) {
		for (subject in 1:length(probmas)) {
			if (row %in% rownames(probmas[[subject]])) {
				countVec[row] <- countVec[row] + counts[[subject]][row];
				nfish[row] <- nfish[row] + 1;
				for (col in colnames(probmas[[subject]])) {
					probMat[row,col] <- probMat[row,col] + probmas[[subject]][row,col];
				}
			}
		}
	}
	
	nfish[which(nfish == 0)] <- 1; # to avoid dividing by zero
	
	if (byTotal) {
		probMat = probMat / length(probmas);
	} else {
		for (col in colnames(probMat)) {
			probMat[,col] <- probMat[,col] / nfish;
		}
	}
	return(list(probMat = probMat, counts = countVec, nfish = nfish));
}

# Compares the transitional probabilities of every pair of behaviors between the first two experimental groups of <data>.
# These values are compared with three different tests: wilcox.test(), t.test(), and bootstrap2independent(). Transitional
# probabilities are only compared for behaviors that occur in at least <minNumLogs> score logs in each group, and where at
# least one animal had a nonzero transitional probability. Graphs output by the bootstrap function are also saved.
#
# TODO check rows & cols add up to 1.
# > apply(probMatsByGroup$nameOfAGroup$probMat,1,sum)
# or
# > sum(probMatsByGroup$nameOfAGroup$probMat) (byTotal true)
.compareTransitionalProbabilities = function(data, byTotal = FALSE, outfilePrefix,
											 tests = list(t.test = t.test, wilcox = wilcox.test, bootstrap = list(func = .bootstrapWrapper)), minNumLogs = 3) {
	data = .filterDataList(data, renameStartStop = TRUE);
	groupwiseLogs = .sepGroups(data);
	
	transProbsByGroup = list();
	for (group in groupwiseLogs$groupNames) {
		groupPMs = lapply(groupwiseLogs$groupData[[group]], function(d) {.getProbabilityMatrix(d$behavior, byTotal=byTotal)});
		transProbsByGroup[[group]] <- .makeTPMatrix(groupPMs, groupwiseLogs$behnames, byTotal);
		write.csv(transProbsByGroup[[group]], file = paste(outfilePrefix, group, "transitionalprobabilities_datadump.csv", sep = "_"));
		groupPMsAndCounts = list(probMats = groupPMs, counts = lapply(groupwiseLogs$groupData[[group]], function(d) {table(d$behavior)}));
		probMat = .combineProbabilityMatrices(groupPMsAndCounts, groupwiseLogs$behnames, byTotal)$probMat;
		# if ((byTotal && sum(probMat) != 1) || (!byTotal && sum(apply(probMat,1,sum) != 1) > 0)) {
			# print(probMat);
			# print(sum(probMat));
			# print(apply(probMat, 1, sum));
			# warning("Something is wrong about this probability matrix - the probabilities don't add up to 1. See Katrina for help.");
		# } TODO TODO FIX THISSSSS
		write.csv(probMat, file = paste(outfilePrefix, group, "transitionalprobabilities_average.csv", sep = "_"));
	}
	return(.runStats(dataByGroup = transProbsByGroup, outfilePrefix = paste(outfilePrefix, "transitionalprobabilities", sep = "_"),
			tests = tests, minNumLogsForComparison = minNumLogs, skipNA = !byTotal));
}


# Helper function for .compareTransitionalProbabilities(). Returns a matrix with a row for each pair of
# behaviors and a column for each subject, giving the transitional probability for that pair of behaviors
# in that subject. If the lead behavior never occurs in the given subject, and byTotal = FALSE, the value is NA.
.makeTPMatrix = function(groupPMs, behnames, byTotal) {
	behcombos <- character();
	for (beh in behnames) {
		behcombos <- c(behcombos, paste(behnames, "->", beh));
	}
	
	TPmat = if (!byTotal) {matrix(nrow = length(behcombos), ncol = length(groupPMs), dimnames = list(behcombos, names(groupPMs)))}
		  else {matrix(data = 0, nrow = length(behcombos), ncol = length(groupPMs), dimnames = list(behcombos, names(groupPMs)))};
	for (behCombo in behcombos) {
		leaderBeh = gsub(" -> .*", "", behCombo);
		followerBeh = gsub(".* -> ", "", behCombo);
		for (subject in names(groupPMs)) {
			if (leaderBeh %in% dimnames(groupPMs[[subject]])[[1]]) {
				if (followerBeh %in% dimnames(groupPMs[[subject]])[[2]]) {
					TPmat[behCombo, subject] = groupPMs[[subject]][leaderBeh, followerBeh];
				} else {
					TPmat[behCombo, subject] = 0;
				}
			}
		}
	}
	return(TPmat);
}

# Source: ethograms_from_scorevideo.R
# Reads in a probability matrix and returns a list containing:
#    1. $probMat, the input probability matrix
#    2. $hMat, a matrix of entropy
#    3. $h, a vector of the raw entropy for each behavioral code
#    4. $h_norm, a vector of the normalized entropy for each behavioral code
#    5. $h_max, a value representing the maximum possible entropy (which was used to normalize
#        the entropy values)
.computeEntropyProbMatrix = function (probMat) {
	hMat = probMat;
	for (row in 1:nrow(probMat)) {
		hMat[row, ] = .computeEntropyOneState(probMat[row, ]);
	}
	num_states = ncol(probMat);
	h_max = num_states * ( (-1/num_states)*log2(1/num_states) );
	h = apply(hMat, 1, sum);
	h_norm = h / h_max;
	return(list(probMat=probMat, hMat=hMat, h=h, h_norm=h_norm, h_max=h_max));
}

# Source: ethograms_from_scorevideo.R
# Helper function for .computeEntropyProbMatrix
.computeEntropyOneState = function (probVec) {
	h = c();
	for (follower in 1:length(probVec)) {
		prob = probVec[follower];
		if (prob == 0) {
			h = c(h, 0);
		} else {
			h = c(h, (-prob)*log2(prob));
		}
	}
	names(h) = names(probVec);
	return(h);
}


# Calls .combineEntropyVecs on each experimental group in <data>, and returns a list of the group-level
# entropy vectors for each group. Notice that the combined entropy is (entropy for 1st subj + entropy for 2nd subj + ...) / (number of subjects),
# NOT the entropy of the combined probability matrix for each group.
.groupLevelEntropy = function(data) {
	data = .filterDataList(data, renameStartStop = TRUE);
	groupwiseLogs = .sepGroups(data);
	# 	return(list(groupNames = groupNames, groupData = groupData, behnames = behnames));
	
	entropyMatsByGroup = list();
	for (group in groupwiseLogs$groupNames) {
		probMats = lapply(groupwiseLogs$groupData[[group]], function(d) {.getProbabilityMatrix(d$behavior, byTotal=FALSE)});
		entropyVecs = lapply(probMats, function(d) {.computeEntropyProbMatrix(d)$h_norm});
		entropyMatsByGroup[[group]] <- .combineEntropyVecs(entropyVecs, groupwiseLogs$behnames);
	}
	return(entropyMatsByGroup);
}

# Helper function for .groupLevelEntropy()
# Takes <entropyVecs>, a list of the entropy vector for each behavior for each subject, and combines them
# by averaging into a unified entropy vector.
.combineEntropyVecs = function(entropyVecs, behnames) {
	#print(entropyVecs);
	entropies = numeric(length = length(behnames));
	names(entropies) <- behnames;
	nfish <- entropies;
	for (beh in names(entropies)) {
		for (subject in 1:length(entropyVecs)) {
			if (beh %in% names(entropyVecs[[subject]])) {
				entropies[beh] = entropies[beh] + entropyVecs[[subject]][beh];
				nfish[beh] <- nfish[beh] + 1;
			}
		}
	}
	
	nfish[which(nfish == 0)] <- 1; # to avoid dividing by zero
	
	return(entropies / nfish);
}

# Compares the entropies of every behavior between the first two experimental groups of <data>. These values
# are compared with three different tests: wilcox.test(), t.test(), and bootstrap2independent(). Entropies are
# only compared for behaviors that occur in at least <minNumLogs> score logs in each group. Graphs output by
# the bootstrap function are also saved.
# TODO combine with .calcBasicStats???
.compareEntropy = function(data, outfilePrefix, tests = list(t.test = t.test, wilcox = wilcox.test,
						   bootstrap = list(func = .bootstrapWrapper)), minNumLogs = 3) {
	data = .filterDataList(data, renameStartStop = TRUE);
	groupwiseLogs = .sepGroups(data);
	# 	return(list(groupNames = groupNames, groupData = groupData, behnames = behnames));
	
	entropiesByGroup = list();
	for (group in groupwiseLogs$groupNames) {
		probMats = lapply(groupwiseLogs$groupData[[group]], function(d) {.getProbabilityMatrix(d$behavior, byTotal=FALSE)});
		entropyVecs = lapply(probMats, function(d) {.computeEntropyProbMatrix(d)$h_norm});
		# groupEMsAndCounts = list(probMats = lapply(entropyLists, function(l) {return(l$hMat)}),
								 # counts = lapply(groupwiseLogs$groupData[[group]], function(d) {table(d$behavior)}));
		entropiesByGroup[[group]] <- .makeEntropyVecMatrix(entropyVecs, groupwiseLogs$behnames);
		write.csv(entropiesByGroup[[group]], file = paste(outfilePrefix, group, "entropydata.csv", sep = "_"));
	}
	return(.runStats(dataByGroup = entropiesByGroup, outfilePrefix = paste(outfilePrefix, "entropy", sep = "_"),
					 tests = tests, minNumLogsForComparison = minNumLogs, skipNA = T));
}


# Helper function for .compareEntropy(). Returns a matrix with a row for each behavior and a column
# for each subject, giving the entropy for that behavior in that subject. If the behavior never occurs
# in the given subject, the value is NA.
.makeEntropyVecMatrix = function(entropyVecs, behnames) {
	entropies = matrix(nrow = length(behnames), ncol = length(entropyVecs), dimnames = list(paste(behnames, ": Entropy"), names(entropyVecs)));
	for (beh in 1:length(behnames)) {
		for (subject in names(entropyVecs)) {
			if (behnames[beh] %in% names(entropyVecs[[subject]])) {
				entropies[beh, subject] = entropyVecs[[subject]][behnames[beh]];
			}
		}
	}
	return(entropies);
}



#####################################################################################################
## MARKOV CHAINS                                                                                   ##
#####################################################################################################

# Source: behavior_syntax.R
# Given a probability matrix and a data vector, writes a dot file to <file> (default is to console)
# Options:
#	minValForLine - only draw lines with probability/frequency greater than minValForLine. Should be between 0 and 1.
#	weird - hacky fix for drawing markov chains weighted by time, where smaller numbers should correspond to thicker lines. Probably don't use this?
#	singleCharLables - puts labels inside the circles that are large enough to hold a single 24-pt character. Default is all labels outside.
#	byTotal - was byTotal on or off when creating the probability matrix? (used in line weighting)#
.buildDotFile = function (probMatrix, originalDataVec, file='', title='untitled', fontsize=24, minValForLine = 0, weird = FALSE, singleCharLabels = FALSE, byTotal = FALSE) {
	# write top line to file
	cat('digraph', title, '\n', '	{\n', file=file);
		
	# # get behavior frequencies
	if (is.numeric(originalDataVec)) {freqs = originalDataVec;}
	else if (is.character(originalDataVec)) {freqs = table(originalDataVec);}
	else {stop("DATA VEC IS WRONG IN .buildDotFile!!!!!!!");}	
	
	#print(freqs);print(probMatrix)
	# # check that behaviors are in same order in freqs and probMatrix
	if (!(sum(rownames(probMatrix)==names(freqs)) == length(freqs))) {stop('NAMES DONT MATCH, GO FIND AUSTIN!!!')}
	
   	toRemove = which(apply(probMatrix, 1, sum) == 0);
	if(length(toRemove)) {
		probMatrix = probMatrix[-toRemove,-toRemove];
		freqs = freqs[-toRemove];
    }

	# # compute proportions of behaviors for relative node size
	for (beh in 1:length(freqs))
	{
#		if (names(freqs)[beh] %in% c("[","]")) {next;}
		prop = freqs[beh] / sum(freqs) * 10; #print(file)
		if (!singleCharLabels || prop < 0.7) { # 0.7 is the magic size for single characters in 24pt font
			cat('		', gsub('[^A-Za-z1-9]', '', names(freqs)[beh]), ' [label="", xlabel="', gsub(' ', '', names(freqs)[beh]),'", width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
		} else {
			cat('		', gsub('[^A-Za-z1-9]', '', names(freqs)[beh]), ' [width=', prop, ', height=', prop, ', fontsize=', fontsize, '];\n', file=file, append=T, sep='');
		}
	}
	
	# loop through probMatrix to get probabilities
	# probMatrix=probMatrix*10;
	for (row in 1:nrow(probMatrix))
	{
		for (col in 1:ncol(probMatrix))
		{
			val = if (byTotal) {(probMatrix[row,col] / sum(probMatrix)) * 100} else if (!is.nan(probMatrix[row,1])) {(probMatrix[row,col] / sum(probMatrix[row,])) * 10} else {-1};
			prob = if(byTotal) val / 100 else val / 10; 
			if (weird) {val = 10 * (1 - (probMatrix[row,col] / max(probMatrix[row,])));} #FIX THIS !!!TODO
			if (prob > minValForLine)
			{
				leader = rownames(probMatrix)[row];
				follower = colnames(probMatrix)[col];#print(paste(leader, follower, sep = ','));
				if (leader == "STOP") {next;}
#				if (follower %in% c("[","]")) {next;}
		
				cat('		', gsub('[^A-Za-z1-9]', '', leader), ' -> ', gsub('[^A-Za-z1-9]', '', follower),
				    ' [label="", style="setlinewidth(', val, ')", arrowsize=1];','\n' ,sep='', file=file, append=T);	
		 	}
		}
	}
	
	# write last line of file
	cat('	}', file=file, append=T);
	
	return(NULL)
}

# Takes a data frame, converts the descriptions to single letters using key (or using
# a smart algorithm to assign letters if no key is provided), and returns a list of the
# data frame and the key used. A key should be a list whose names are the behavior names
# and whose entries are single-letter codes.
.descriptionsToLetters = function (data, key = NULL) {
	frequencies = table(data$behavior);
	if (is.null(key)){key = .assignLetters(frequencies);}
	else {.validateKey(key, names(frequencies))}
	for (i in 1:dim(data)[1]) {
		data$behavior[i] <- key[[data$behavior[i]]];
	}
	descTable = data.frame(code = as.character(key), description = names(key));
	#dimnames(descTable)[[1]] <- 1:dim(descTable)[1];
	return(list(data, descTable));
}

# Helper function for .descriptionsToLetters. Checks that <key> is valid.
.validateKey = function (key, behaviors) {
	if (!is.list(key)) stop("Error: Key provided to .descriptionsToLetters is not a list.\nKey should be a list of single-letter codes, with the names of key set to the behavior names.");
	if (sum(!(behaviors %in% names(key))) != 0) stop("Error: Key provided to .descriptionsToLetters does not contain a code for all behaviors in the data provided.");
	asvector = as.character(key);
	if(sum(nchar(asvector) != 1) != 0) {warning("Not all codes in key are single characters");}
}

# Produces a key for .descriptionsToLetters if no key is provided.
.assignLetters = function (freqtable) {
	behList = names(freqtable[order(freqtable, decreasing = TRUE)]);
	key = list();
	orphans = character();
	for (behaviorName in behList) {
		preferredCodes = toupper(substr(unlist(strsplit(behaviorName, " ")), 1,1));
		for (code in preferredCodes) {
			if ( !(code %in% as.character(key)) ) {
				key[[behaviorName]] <- code;
				break;
			}
		}
		if (! behaviorName %in% names(key)) {orphans <- c(orphans, behaviorName);}
	}
	
	for (behaviorName in orphans) {
		for (code in c(toupper(letters), letters)) {
			if ( !(code %in% as.character(key)) ) {
				key[[behaviorName]] <- code;
				break;
			}
		}
		if (! behaviorName %in% names(key)) {stop("More than 52 behaviors are present. .assignLetters() needs modification.")}
	}
	return(key);
}


# Makes a transition probability matrix for <data>, either <byTotal> or not, then draws a markov chain in file <filename> with
# line cutoff <minValForLine>. Additional arguments (...) are passed to a call to .filterData before any calculations take place.
# OPTIONS:
# minValForLine - minimum probability*10 value to justify drawing an arrow
# k - (ARCHAIC) number of spots that are ok to count for transitional probabilities. weighted makes behaviors closer to leader count more.
.makeDotPlot = function (data, filename, byTotal = FALSE, minValForLine = 0, singleLetterLabels = FALSE, ...) {
	data = .filterData(data, ...);
	
	cleanerDataForPlot <- .renameStartStop(data);
	behaviorVec <- cleanerDataForPlot$behavior;
	probMatrix <- .getProbabilityMatrix(behaviorVec, byTotal=byTotal);
	.buildDotFile(probMatrix, behaviorVec, file = filename, minValForLine = minValForLine);
	return(list(probMatrix = probMatrix, entropyData = .computeEntropyProbMatrix(probMatrix)));
	
	# cleanerDataForPlot <- .descriptionsToLetters(.renameStartStop(data));
	# behaviorVec <- cleanerDataForPlot[[1]]$behavior;
	# probMatrix <- .getProbabilityMatrix(behaviorVec, byTotal=byTotal);
	# .buildDotFile(probMatrix, behaviorVec, file = filename, minValForLine = minValForLine);
	# return(list(probMatrix = probMatrix, codes = cleanerDataForPlot[[2]], entropyData = .computeEntropyProbMatrix(probMatrix)));
}

# Calls .makeDotPlot on every data frame in <listOfData>.
# TODO make a list of the lists returned by .makeDotPlots. so entropy, probmatrix, codes, etc are not lost.
.makeDotPlots = function (listOfData, fileprefix, filename_ext = ".dot", ...) {
	filenames = paste(fileprefix, gsub('/', '.', names(listOfData)), filename_ext, sep = '_')
	for (i in 1:length(filenames)) {
		.makeDotPlot(listOfData[[i]], filenames[i], ...)
	}
}

# Makes markov chains from the data format output by .groupLevelProbMats().
.makeDotPlotsFromProbMas = function(probMatDat, fileprefix, ...) {
	for (i in 1:length(probMatDat)) {
		.buildDotFile(probMatDat[[i]]$probMat, probMatDat[[i]]$counts, file = paste(fileprefix, names(probMatDat)[[i]], ".dot", sep = ""), ...);
	}
}

# Makes group markov chains. Calls .filterDataList with the ... arguments, then calls .groupLevelProbMats, and finally calls .makeDotPlotsFromProbMas.
.makeGroupDotPlots = function(data, fileprefix, byTotal = FALSE, minValForLine = 0, singleLetterLabels = FALSE, ...) {
	data = .filterDataList(data, ...);
	
	cleanerDataForPlot <- .filterDataList(data, renameStartStop=T);
	probMatData <- .groupLevelProbMats(cleanerDataForPlot, byTotal=byTotal);
	.makeDotPlotsFromProbMas(probMatData, fileprefix, byTotal=byTotal, minValForLine=minValForLine, singleCharLabels=singleLetterLabels);
}


#####################################################################################################
## COLOR KEYS                                                                                      ##
#####################################################################################################

# Validates the behaviorsToPlotAndColors provided to .behavioralDensityGraph() or .behavioralDensityGraphs()
.validateColorKey = function(behcolors, validBehNames = NULL) {
	if (dim(behcolors)[2] != 2) stop("behaviorsToPlotAndColors must have 2 columns!");
	if (sum(!(behcolors[,2] %in% colors())) != 0) {
		stop(paste('Invalid color in behaviorsToPlotAndColors: "', behcolors[!(behcolors[,2] %in% colors()), 2], '"\n', sep = ""));
	}
	if (!is.null(validBehNames) && sum(!(behcolors[,1] %in% validBehNames)) != 0) {
		stop(paste('Behavior in behaviorsToPlotAndColors: "', behcolors[!(behcolors[,1] %in% validBehNames), 1], '" not found in any score log\n', sep = ""));
	}
}

# Plots a color key to make it human-readable
.plotColorLegend = function(colorkey) {
	behaviors = colorkey[,1];
	colors = colorkey[,2];
	par(mai = c(1, 3, 1, 1))
	plot(c(0, 1), c(0, length(behaviors)), frame.plot=F, axes=F, xlab = '', ylab='', xlim = c(0, 1), ylim = c(0.5, length(behaviors) + 0.5), col = "white", main = "Color Key");
	for (i in 1:length(behaviors)) {
		rect(xleft = 0, xright = 1, border = NA,
			 col = colors[i], ybottom = length(behaviors) - i + .5, ytop = length(behaviors) - i+1.5);
	}
	axis(2, at=1:length(behaviors), labels=behaviors[length(behaviors):1], tick=F, las=2);
}

# Prints out the color key in a nice format to the console, with each row preceded by <whitespace>.
# Also plots it.
.printColorKey = function(colorKey, whitespace = "") {
	paddingLengths = 1 + max(nchar(colorKey)) - nchar(colorKey);
	for (i in 1:length(colorKey[,1])) {
		cat(whitespace, '"', colorKey[i,1], '":', rep_len(' ', paddingLengths[i]), '"', colorKey[i, 2], '"\n', sep = "");
	}
	.plotColorLegend(colorKey);
}

# Prompts the user to enter a color for <beh> and reprompts until they enter either a valid color or "none".
.getColor = function(beh) {
	prompt = paste("  What color should \"", beh, '" be? (Enter "none" to not plot this behavior)\n  > ', sep = "");
	userInput = gsub('^["\']','', gsub('["\']$','', gsub(' ', '', readline(prompt))));
	userInput = .autocomplete(userInput, c("none", colors()));
	while (!(userInput %in% c("none", colors()))) {
		prompt = paste("  Please enter a valid color name or \"none\", or press ESC to quit. What color should \"", beh, "\" be?\n  > ", sep = "");
		userInput = gsub('^["\']','', gsub('["\']$','', gsub(' ', '', readline(prompt))));
		userInput = .autocomplete(userInput, c("none", colors()));
	}
	return(userInput);
}

# Allows rows to be added to <colorkey>, removed from <colorkey>, or have their colors changed.
.editColorKey = function(colorkey, validBehNames) {
	prompt = "Enter a behavior name to edit its color, \"p\" to print the current key, or \"q\" to quit color editor.\n  > ";
	userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
	userInput = .autocomplete(userInput, c("q", "p", validBehNames));
	while (userInput != "q") {
		if (userInput %in% validBehNames) {
			color = .getColor(userInput);
			if (color == "none" && userInput %in% colorkey[,1]) {
				colorkey <- colorkey[-which(colorkey[,1] == userInput),];
			} else if (userInput %in% colorkey[,1]) {
				colorkey[which(colorkey[,1] == userInput), 2] <- color;
			} else {
				colorkey = rbind(colorkey, c(userInput, color));
			}
			.plotColorLegend(colorkey);
		} else if (userInput == "p") {
			.printColorKey(colorkey, "  ");
		}
		userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
		userInput = .autocomplete(userInput, c("q", "p", colorkey[,1]))
	}
	return(colorkey);
}

# Allows the user to reorder <colorkey> by typing in the behaviors in a new order.
.reorderColorKey = function(colorkey) {
	prompt = "  Which behavior do you want in the very back? (type \"l\" to list remaining behaviors, or \"q\" to cancel reordering)\n  > ";
	userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
	userInput = .autocomplete(userInput, c("l", "q", colorkey[,1]));
	prompt = "  Which behavior do you want next? (type \"l\" to list remaining behaviors, or \"q\" to cancel reordering)\n  > ";
	newcolorkey = NULL;
	behsleft = colorkey[,1];
	while (length(behsleft) > 0) {
		if (userInput %in% behsleft) {
			newcolorkey = rbind(newcolorkey, colorkey[which(colorkey[,1] == userInput),]);
			behsleft = behsleft[-which(behsleft == userInput)];
		} else if (userInput == "l") {
			cat('  "');
			cat(behsleft, sep = '" "');
			cat('"\n');
		} else if (userInput == "q") {
			return(colorkey);
		} else {
			if (userInput %in% newcolorkey[,1]) {
				cat("  Behavior \"", userInput, "\" was already added. ", sep = "");
				if (.getYesOrNo("  Would you like to move it to the front? ")) {
					newcolorkey = rbind(newcolorkey[-which(newcolorkey[,1] == userInput),], newcolorkey[which(newcolorkey[,1] == userInput),]);
				}
			} else {
				cat("  Not a valid behavior.\n  Behaviors left:\n    \"");
				cat(behsleft, sep = '" "');
				cat('"\n');
			}
		}
		
		if (length(behsleft > 0)) {
			userInput = gsub('^["\']','', gsub('["\']$','', readline(prompt)));
			userInput = .autocomplete(userInput, c("l", "q", colorkey[,1]));
		}
	}
	return(newcolorkey);
}

# A user-friendly color key constructor.
.buildColorKey = function(validBehNames) {
	cat("Building color key...\nI found these behaviors:\n\"");
	cat(validBehNames, sep = '" "');
	cat('"\n');
	
	colorkey = NULL;
	for (beh in validBehNames) {
		color = .getColor(beh);
		if (color != "none") {
			colorkey = rbind(colorkey, c(beh, color));
		}
	}
	
	cat("Color key created: \n");
	.printColorKey(colorkey, "  ");
	while (!(.getYesOrNo("Are these colors okay? "))) {
		colorkey = .editColorKey(colorkey, validBehNames);
		cat("Color key created: \n");
		.printColorKey(colorkey, "  ");
	}
	
	cat("In graphs, behaviors will be plotted in the order shown in this key.\nThat is, the behaviors at the bottom of the key will be plotted in front of those at the top.\n")
	while (!(.getYesOrNo("Is this order okay? "))) {
		colorkey = .reorderColorKey(colorkey);
		cat("Color key created: \n");
		.printColorKey(colorkey, "  ");
	}

	return(colorkey);
}


#####################################################################################################
## BEHAVIORAL DENSITY PLOTS                                                                        ##
#####################################################################################################

# Returns a vector giving the distances between occurances of <centerBeh> and <varBeh> in data frame <data>.
# If <noRepCenterBeh> is false, the vector has length count(centerBeh)*count(varBeh), and it has distances between every occurance of each behavior.
# If <noRepCenterBeh> is true, the vector has length roughly 2*count(varBeh), as it only measures distances between a varBeh and its two neighboring centerBehs.
# Optionally you can provide a buffer, specifying not to count centerBehs that occur too close to the beginning or end of data. This buffer
# can be given in time (seconds) or behavior counts.
# This function returns a list containing:
#     1. <timeDists>, a vector of the distances in times
#     2. <behDists>, a vector of the distances in behaviors
#     3. <centerCount>, the number of occurances of behavior centerBeh
#     4. <varCount>, the number of occurances of behavior varBeh
.getAllIntervals = function (data, centerBeh, varBeh, startBuffer = 0, endBuffer = 0, bufferInTimes = TRUE, noRepCenterBeh = TRUE) {
	varBehLocs = which(data$behavior == varBeh);
	varBehTimes = data$time[varBehLocs];
	centerBehLocs = which(data$behavior == centerBeh);
	if (bufferInTimes) {
		max_time = max(data$time);
		centerBehTimes = data$time[centerBehLocs];
		centerBehLocs = centerBehLocs[centerBehTimes > startBuffer & centerBehTimes <= max_time - endBuffer];
	} else {
		centerBehLocs = centerBehLocs[as.numeric(centerBehLocs) > startBuffer & as.numeric(centerBehLocs) <= length(data$behavior - endBuffer)];
	}
	centerBehTimes = data$time[centerBehLocs];
	
	
	behDists = numeric();
	timeDists = numeric();
	
	for (i in 1:length(centerBehLocs)) {
		if (noRepCenterBeh) {
			if (i == 1) {
				behDists = c(behDists, varBehLocs[varBehLocs <= centerBehLocs[i+1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes <= centerBehTimes[i+1]] - centerBehTimes[i]);				
			} else if (i == length(centerBehLocs)) {
				behDists = c(behDists, varBehLocs[varBehLocs >= centerBehLocs[i-1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes >= centerBehTimes[i-1]] - centerBehTimes[i]);				
			} else {
				behDists = c(behDists, varBehLocs[varBehLocs >= centerBehLocs[i-1] & varBehLocs <= centerBehLocs[i+1]] - centerBehLocs[i]);
				timeDists = c(timeDists, varBehTimes[varBehTimes >= centerBehTimes[i-1] & varBehTimes <= centerBehTimes[i+1]] - centerBehTimes[i]);
			}
		} else {
			behDists = c(behDists, varBehLocs - centerBehLocs[i]);
			timeDists = c(timeDists, varBehTimes - centerBehTimes[i]);
		}
	}
	
	return(list(behDists=behDists, timeDists=timeDists, centerCount = length(centerBehLocs), varCount = length(varBehLocs)));
}

# Calls .getAllIntervals on every data frame in the list <data>, and combines the results into one
# list. behDists and timeDists contain all the distances from all the fish, and varCount and centerCount
# contain the total counts across all fish.
.getIntervalsAcrossFish = function (data, ...) {
	alldata <- list(behDists = numeric(), timeDists = numeric(), centerCount = 0, varCount = 0);
	for (i in 1:length(data)) {
		tmp <- .getAllIntervals(data[[i]], ...);
		alldata$behDists <- c(alldata$behDists, tmp$behDists);
		alldata$timeDists <- c(alldata$timeDists, tmp$timeDists);
		alldata$centerCount <- alldata$centerCount + tmp$centerCount;
		alldata$varCount <- alldata$varCount + tmp$varCount;
	}
	return(alldata);
}


# Makes a behavioral density plot of the behaviors around <centerBeh> for either a single fish (multifish = FALSE) or every fish in
# a list (multifish = TRUE). Each behavior to be plotted must have a row in <behaviorsToPlotAndColors>: its name in column 1, and the color
# it should be in column 2. It is the caller's responsibility to make this key and remember which color corresponds to which behavior.
# If a <filename> is given, the plot is saved to that file. <lim> represents the limits of the time axis in seconds; the default (15) gives
# fifteen seconds before and after <centerBeh>. If <noRepCenterBeh> is TRUE, the plotted behaviors are only counted in each direction of a
# <centerBeh> if another <centerBeh> has not yet been reached. If it is FALSE, the plotted behaviors are counted until the beginning and end
# of the assay. This parameter will make no difference for behaviors that occur wide apart, but it may make a difference for behaviors
# that are spaced closely together. <timesPerBin> controls the amount of time that is binned together (default 0.5 seconds) and ymax
# gives the maximum value of the y-axis.
# weightingStyle refers to the way the y axis is supposed to be normalized. If it is "singlebeh", the y-values represent the fraction of
# the total occurrances of the plotted behavior that occur in a given time bin. If it is "allbeh". the y-values represent the fraction of
# all behaviors that are (1) the plotted behavior and (2) in the given time bin. If it is "rawcounts", the y-value is just the count of
# the plotted behavior that occurs in that timebin. # TODO normalize by centering beh
# Note that if you call this function directly, there is no check on the names of behaviors in behaviorsToPlotAndColors. This enables
# you to do things like give the same color key for each group even if a behavior in it is never performed in a given group, but it
# also enables you to make errors. Be careful!
# TODO hist & normalize for each fish separately, THEN combine
# TODO sliding window
# TODO write function that makes generation of <behaviorsToPlotAndColors> easier
# TODO write a function that automatically generates a nice color key given <behaviorsToPlotAndColors>
# TODO get the data out in a nice way (aka not a chart)
.behavioralDensityGraph = function(data, behaviorsToPlotAndColors, centerBeh, filename = NULL, lim = 15, noRepCenterBeh = TRUE, multifish = FALSE,
								   timesPerBin = 0.5, ymax = 1, weightingStyle = "singlebeh") {
	.validateColorKey(behaviorsToPlotAndColors);
	data <- if (multifish) .filterDataList(data, renameStartStop = TRUE) else .renameStartStop(data);
	if (!is.null(filename)) jpeg(filename = filename, width = 15, height = 5, units = "in", quality = 100, res = 150, type = "quartz");
	histBreaks = ((-ceiling(lim/timesPerBin) - 1):(ceiling(lim/timesPerBin)) + 0.5) * timesPerBin;
	for (i in 1:(dim(behaviorsToPlotAndColors)[1])) {
		x <- if (multifish) {.getIntervalsAcrossFish(data, centerBeh = centerBeh, varBeh = behaviorsToPlotAndColors[i, 1], noRepCenterBeh = noRepCenterBeh);}
			 else {.getAllIntervals(data, centerBeh, behaviorsToPlotAndColors[i, 1], noRepCenterBeh = noRepCenterBeh);}
		h <- hist(x$timeDists[x$timeDists > -lim & x$timeDists < lim], breaks = histBreaks, plot = FALSE);
		if (weightingStyle == "singlebeh") {
			count <- x$varCount;
			plot(x = h$mids, y = h$counts / count, type = "l", col = behaviorsToPlotAndColors[i, 2], xlim = c(-lim, lim), ylim = c(0, ymax), xlab = "", ylab = "");
		} else if (weightingStyle == "allbeh") {
			count <- if (multifish) {sum(unlist(lapply(data, function (d) {length(d$behavior)})));}
			         else {length(d$behavior);}
			plot(x = h$mids, y = h$counts / count, type = "l", col = behaviorsToPlotAndColors[i, 2], xlim = c(-lim, lim), ylim = c(0, ymax), xlab = "", ylab = "");
		} else if (weightingStyle == "rawcounts") {
			plot(x = h$mids, y = h$counts, type = "l", col = behaviorsToPlotAndColors[i, 2], xlim = c(-lim, lim), ylim = c(0, ymax), xlab = "", ylab = "");
		} else {
			stop("ERROR IN .behavioralDensityGraph(): Invalid weightingStyle.");
		}
		#print(paste("mids length =", length(h$mids), ", counts length =", length(h$counts), ", count =", count, ",extra length =", length(h$counts / count)))
		par(new = TRUE);
	}
	title(main = centerBeh, xlab = paste("Time after", centerBeh, "(seconds)"),
			ylab = if(weightingStyle=="singlebeh") "Fraction of individual behavior" else if(weightingStyle=="allbeh") "Fraction of all behaviors" else "Count");
	par(new = FALSE);
	dev.off();
}

# Draws behavioral density graphs for each experimental group in <data>. If no targetBehs are provided, a behavioral density plot is made
# centered at each behavior in behaviorsToPlotAndColors; if targetBehs are provided, graphs are only made for the behaviors in it.
# For example, calling with targetBehs = c("Lead", "Quiver") would make a graph centered at "Lead" for each group and a graph centered
# at "Quiver" for each group.
# For more information, look at the documentation for .behavioralDensityGraph().
.behavioralDensityGraphs = function(data, behaviorsToPlotAndColors, filePref, targetBehs = NULL, ...) {
	data = .filterDataList(data, renameStartStop = TRUE);
	groupwiseLogs = .sepGroups(data);

	if (is.null(targetBehs)) {behnames = behaviorsToPlotAndColors[,1];}
	else {
		for (beh in targetBehs) {
			if (!(beh %in% groupwiseLogs$behnames)) stop(paste("Error: Behavior", beh, "is not in any score log."))
		}
		behnames = targetBehs;
	}
	
	.validateColorKey(behaviorsToPlotAndColors, groupwiseLogs$behnames);
	
	for (beh in behnames) {
		cat("Plotting behavior \"", beh, '"...\n', sep = "");
		par(mfrow = c(length(groupwiseLogs$groupNames), 1));
		for (i in 1:length(groupwiseLogs$groupNames)) {
			.behavioralDensityGraph(groupwiseLogs$groupData[[i]], behaviorsToPlotAndColors, centerBeh = beh,
								    filename = paste(filePref, "_", beh, "_", groupwiseLogs$groupNames[i], ".jpeg", sep = ""), multifish = TRUE, ...);
		}
#		readline("Press ENTER when you are ready to continue:");
	}
}

#####################################################################################################
## RASTER PLOTS                                                                                    ##
#####################################################################################################

# Calls .makeRasterPlot on every data frame in <data> and saves the results as jpegs
#   to a new folder called "raster_plots" in <outfilePrefix>.
.makeRasters = function(data, outfilePrefix) {
	outdir = paste(outfilePrefix, "raster_plots/", sep = "");
	dir.create(outdir);
	for (i in 1:length(data)) {
		print(names(data[i]));
		filename = paste(outdir, gsub("\\.txt", ".jpeg", gsub("(.*/)*", "", names(data)[i])), sep = "");
		jpeg(filename = filename, width = 7.5, height = 10, units = "in",
		     quality = 100, res = 300, type = "quartz");
		.makeRasterPlot(data[[i]]);
		dev.off();
	}
}

# Source: something from Austin
# Makes a raster plot for a single fish with each behavior on a row by itself. A little bit jenky
# in terms of how the tick marks are drawn. Also makes boxplots that are presumably meaningful
# in some way, and returns something.
.makeRasterPlot = function (dataFrame, plots=T, ...) {
	data = dataFrame$behavior;
	# if (is.null(.checkInputDataVecOK(data))) {
		# return(NULL);
	# } 
	codes = names(table(dataFrame$behavior));
	frames = as.numeric(dataFrame$time);
	counts = table(data);
	num_beh = length(counts);
	diffs = diff(frames);
	
	avg_diffs = counts;
	for (beh in 1:num_beh) {
		avg_diffs[beh] = mean(diffs[data==names(avg_diffs)[beh]], na.rm=T);
	}
	
	if (plots) {
		par(mfrow=c(2, 1));
		par(oma=c(0,5,0,0));
		plot(frames, as.numeric(as.factor(data)), frame.plot=F, 
			 axes=F, xlab='time (seconds)', ylab='', 
			 col='blue', cex=3, pch=3, 
			 ...);
		axis(2, at=1:num_beh, labels=codes, tick=F, las=2);
		axis(1, yaxp=c(0, max(frames), 10), col='white', col.ticks='black');
		for (i in 1:num_beh) { 
			abline(h=i, col='darkgrey');
		}
		
		boxplot(frames ~ as.factor(data), frame.plot=F,
				col='grey', notch=F, width=table(data)/length(data), 
				horizontal=T, names=codes, las=1
				);
	}
	probMat = .getProbabilityMatrix(data);
	return(list(beh_counts=counts, 
				frame_diffs=diffs, 
				frame_diffsAvg=avg_diffs, 
				frames_per_beh=max(frames)/sum(counts),
				ethogram=.computeEntropyProbMatrix(probMat),
				data=dataFrame
				)
			);
}

# Returns a vector of the names of behaviors that are durational in <data>, a list of logs.
.startStopBehs = function(data) {
	behaviors <- names(.findDupBehaviors(data));
	ssBehs = character();
	for (beh in behaviors) {
		if (sum(!is.na(unlist(lapply(data, function(d){d$type[d$behavior == beh][1]})))) > 0) ssBehs <- c(ssBehs, beh);
	}
	return(ssBehs);
}

# Separates the logs in <data> by group and makes a multicolor raster plot of each one. The plots for
# each group are saved to files starting with <outfilePrefix>. Parameters in the ... are passed on to
# .makeMulticolorRasterPlot().
#
# <behaviorsToPlotAndColors> and <durationalBehs> work the same as in .makeMulticolorRasterPlot(). If no
# <behaviorsToPlotAndColors> is provided, the user is prompted to create a color key.
#
# The user can optionally provide a <sortAttribute> to order logs by an attribute (for example, GSI,
# time spent in pot, number of spawns, latency to quiver 5 times, etc.) If this is provided, the parameters
# sort.na.last and sort.decreasing are passed through to the call to order(). Additionally, the value
# will be displayed on the y-axis along with sort.name, which should be something like "GSI", "Time in pot",
# "Spawning count", "Quiver latency", etc.
.makeMulticolorRasterPlots = function (data, outfilePrefix, behaviorsToPlotAndColors = NULL, durationalBehs = NA,
										sortAttribute = NULL, sort.na.last = T, sort.decreasing = F, sort.name = "", ...) {
	if (!is.null(sortAttribute)) {
		names(data) <- paste(names(data), "                               ", sort.name, ":", as.character(sortAttribute));
		data <- .sortByAttribute(data, sortAttribute, na.last = sort.na.last, decreasing = sort.decreasing);
	}
	# rezero if desired here TODO
	
	groupwiseLogs = .sepGroups(data);
	
	if (!is.null(behaviorsToPlotAndColors)) {
		.validateColorKey(behaviorsToPlotAndColors, groupwiseLogs$behnames);
	} else {
		cat("No color key provided to .makeMulticolorRasterPlots. Follow the prompts to create a key, or press ESC at any time to exit.\n")
		behaviorsToPlotAndColors = .buildColorKey(groupwiseLogs$behnames);
	}
	
	startStopBehs = NULL;
	for (dataset in groupwiseLogs$groupData) {
		startStopBehs = c(startStopBehs, .startStopBehs(dataset));
	}
	
	if (is.na(durationalBehs[1])) {
		durationalBehs = names(table(startStopBehs));
	} else {
		if (sum(!(durationalBehs %in% names(table(startStopBehs)))) > 0) {
			stop(paste('Behavior "', durationalBehs[which(!(durationalBehs %in% names(table(startStopBehs))))][1], 
						'" in durationalBehs was not scored as a durational behavior.', sep = ""));
		}
	}
	
	# print(durationalBehs);
	
	for (i in 1:length(groupwiseLogs$groupNames)) {
		.makeMulticolorRasterPlot(groupwiseLogs$groupData[[i]], behaviorsToPlotAndColors,
									filename = paste(outfilePrefix, "_rasterplot_", groupwiseLogs$groupNames[i], ".jpeg", sep = ''),
									durationalBehs = durationalBehs, ...)
	}
}


# Makes a raster plot using color key <behaviorsToPlotAndColors> with one row for each log
# in <data>. If <filename> is not null, the result is saved as a jpg named <filename> that
# has size <widthInInches> x <heightInInches>. Only behaviors in the color key are plotted - 
# to suppress plotting of a given behavior, just exclude it from the color key.
#
# If <durationalBehs> is provided, the behaviors in <durationalBehs> are plotted as bars
# from start to stop instead of ticks. Durational behaviors not in that list but in
# <behaviorsToPlotAndColors> will get a tick for start and a tick for stop. If not provided,
# all behaviors that are durational in the logs will be plotted as bars.
#
# If <staggerSubjects> is set to TRUE, male behaviors are plotted above the line and female
# behaviors below it. This could be rewritten to make it possible to plot an arbitrary list
# of behaviors above/below the line, but right now it is hard-coded in as "male" and "female".
#
# <wiggle> controls the amount of space the durational bars are allowed to take up. It
# should always be a value between 0 and 0.5. The default .2 gives the durational bars 40%
# of the height of the ticks. <defaultDur> is the width, in seconds, of the ticks. If some
# ticks are very light or invisible, try increasing <defaultDur>.
#
# TODO no error checking on durationalBehs. Will probably crash if you put a non-durational beh in
#   there. Also w/ ssBehs, what if it is ss in some-but-not-all logs? Add code to catch this case,
#   draw a line, throw a warning.
# TODO test with 1 durational beh, 1 non-durational beh, 0 of each.
.makeMulticolorRasterPlot = function (data, behaviorsToPlotAndColors, filename = NULL, wiggle = .2, defaultDur = 1,
									  durationalBehs = NA, staggerSubjects = F, widthInInches = 12, heightInInches = 12) {
	if (!is.null(filename)) jpeg(filename = filename, width = widthInInches, height = heightInInches, units = "in", quality = 100, res = 300, type = "quartz");
	
	subjects = names(data);
	num_subj = length(subjects); # TODO change <---- why does this need changing? UNCLEAR.
	maxtime = max(unlist(lapply(data, function(d){max(d$time)})));
	mintime = min(c(0, unlist(lapply(data, function(d){min(d$time)}))));
	ssBehs = if(is.na(durationalBehs[1])) .startStopBehs(data) else durationalBehs;
	
	# print(ssBehs);
	
	.validateColorKey(behaviorsToPlotAndColors);
	singleBehsAndColors = rbind(NULL, behaviorsToPlotAndColors[!(behaviorsToPlotAndColors[,1] %in% ssBehs),]);
	ssBehsAndColors = rbind(NULL, behaviorsToPlotAndColors[behaviorsToPlotAndColors[,1] %in% ssBehs,]);
	# print(singleBehsAndColors);
	# print(ssBehsAndColors);
	
	par(oma=c(0,5,0,0));
	plot(0, 0, frame.plot=F, axes=F, xlab = '', ylab='', xlim = c(mintime, maxtime), ylim = c(0, num_subj + 1), col = "white");

	for (n in 1:num_subj) {
		dataFrame = data[[n]];
		temp_behcolors = singleBehsAndColors[singleBehsAndColors[,1] %in% dataFrame$behavior,];
		for (i in 1:length(temp_behcolors[,1])) {
			beh = temp_behcolors[i,1];
			# print(beh);
			times = dataFrame$time[dataFrame$behavior == beh];
			ybottom = n - 0.5;
			ytop = n + 0.5;
			if (staggerSubjects && dataFrame$subject[dataFrame$behavior == beh][1] %in% c("male", "female")) {
				if (dataFrame$subject[dataFrame$behavior == beh][1] == "male") {
					ybottom = n;
					ytop = n + .45;
				} else {
					ytop = n;
					ybottom = n - .45;
				}
			}
			rect(xleft = times, xright = times + defaultDur,
				 ybottom = ybottom, ytop = ytop,
				 col = temp_behcolors[i,2], border = NA);

			par(new = TRUE);
		}
	}
	par(new = FALSE); 

	title(xlab = "time (seconds)");
	axis(2, at=1:num_subj, labels=subjects, tick=F, las=2);
	axis(1, yaxp=c(mintime, maxtime, 10), col='white', col.ticks='black');
	
	nSsBehs = length(ssBehs);
	durBehBounds = data.frame(bottomBound = ((0:(nSsBehs - 1)) * 2 * wiggle / nSsBehs) - wiggle, topBound = ((1:nSsBehs) * 2 * wiggle / nSsBehs) - wiggle);
	# print(durBehBounds);
	dimnames(durBehBounds)[[1]] <- ssBehsAndColors[,1];
	for (n in 1:num_subj) { 
		abline(h=n, col='black');
		dataFrame = data[[n]];
		temp_behcolors = data.frame(behs = ssBehsAndColors[ssBehsAndColors[,1] %in% dataFrame$behavior,1], colors = ssBehsAndColors[ssBehsAndColors[,1] %in% dataFrame$behavior,2]);
#		print(temp_behcolors);
		for (i in 1:length(temp_behcolors[,1])) {
			beh = temp_behcolors[i,1];
			# print(beh);
			occurrences = data.frame(start = dataFrame$time[dataFrame$behavior == beh & dataFrame$type == "start"],
									 duration = as.numeric(dataFrame$duration[dataFrame$behavior == beh & dataFrame$type == "start"]));
			rect(xleft = occurrences$start, xright = occurrences$start + occurrences$duration,
				 ybottom = n + durBehBounds[beh,]$bottomBound, ytop =  n + durBehBounds[beh,]$topBound,
				 col = temp_behcolors[i,2], border = NA);
		}
	}
	abline(v = 0, col='black');
	
	if (!is.null(filename)) dev.off();
}






###########
# Other Functions
###########






# Reads in a vector giving a sequence of behaviors and returns a matrix giving the probability that the
# behavior in a given column follows the behavior in a fiven row within k behaviors.
.getProbabilityMatrixK = function (data, removeZeroCol=F, k = 2, weighted = FALSE) {
	beh = names(table(data));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeTransitionProbabilityK(data=data, leader=leader, follower=follower, k=k, weighted = weighted);
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$probability;
		}
	}
	if (removeZeroCol) {
		colSums = apply(probMat, 2, sum);
		if (any(colSums==0)) {
			probMat = probMat[-which(colSums==0), -which(colSums==0)];
		}
	}
	return(probMat);
}

# Helper function for .getProbabilityMatrixK
.computeTransitionProbabilityK = function (data, leader, follower, k, weighted) {
	count = 0;
	termination = 0;
	for (i in 1:length(data)) {
		if (data[i] == leader) {
			if (i == length(data)) {
				termination = 1;
			} else {
				increment = 1;
				while (increment <= k && i + increment <= length(data)) {
					if (data[i + increment] == follower) {
						count = if(weighted) {count + (k + 1) - increment} else {count + 1}
					}
					increment = increment + 1;
				}
			}
		}
	}
	total_leader = sum(data == leader);
	prob = count / total_leader;
	return(list(probability=prob, termination=termination, count_transitions=count, count_leader=total_leader));
}

# Reads in a data frame giving a sequence of behaviors with time numbers and returns a list containing matrices:
#   1. timeDistances, the average amount of time between two behaviors,
#   2. behDistances, the average number of behaviors between two behaviors (minimum 1)
#   3. probFollowed, the probability that the leading behavior was followed by the follower at some point. 
# Usage: .getIntervalMatrix(.renameStartStop(dataFrame))
.getIntervalMatrix = function (data) {
	beh = names(table(data$behavior));
	timeMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	behMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	probMat = matrix(nrow=length(beh), ncol=length(beh), dimnames=list(beh, beh));
	for (leader in rownames(probMat)) {
		for (follower in colnames(probMat)) {
			tmp = .computeAverageInterval(data=data, leader=leader, follower=follower);
			timeMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$average_times;
			behMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] = tmp$average_behaviors;
			probMat[match(leader, rownames(probMat)), match(follower, colnames(probMat))] =
				(tmp$count_with / (tmp$count_with + tmp$count_without + tmp$count_end));
		}
	}
	return(list(timeDistances = timeMat, behDistances = behMat, probFollowed = probMat));
}


# Helper function for .getIntervalMatrix
# Computes average time interval, count interval, and returns counts of different situations.
# count_without is the number of times leader was followed by leader before it was followed by follower
# count_end is 1 if there is no follower after the last leader, and 0 otherwise
# count_with is the number of times leader was followed by follower before the end of behavior and before
#       another leader occurred
.computeAverageInterval = function (data, leader, follower) {
	count = 0;
	countWithout = 0; 
	totaltimes = 0;
	totalBehaviors = 0;
	for (i in 1:dim(data)[1]) {
		if (data$behavior[i] == leader) {
			lastStartTime = as.numeric(data$time[i]);
			j = i+1;
			while (j <= dim(data)[1] && data$behavior[j] != follower && data$behavior[j] != leader) {j = j + 1;}
			if (j > dim(data)[1]) {
				countEnd = countEnd + 1;
			} else if (data$behavior[j] == follower) {
				count = count + 1;
				totalTimes = totalTimes + (data$time[j] - lastStartTime);
				totalBehaviors = totalBehaviors + (j-i);
			} else if (data$behavior[j] == leader) {
				countWithout = countWithout + 1;
			} else {
				stop('Something is wrong in ".computeAverageInterval()". GO FIND KATRINA!!!')
			}
		}
	}
	if (count > 0) {
		avgF = totalTimes / count;
		avgB = totalBehaviors / count;
	} else {
		avgF = NA;
		avgB = NA;
	}
	return(list(average_times=avgF, average_behaviors = avgB, count_with = count, count_without=countWithout, count_end=countEnd));
}


# TODO figure out how to preserve time numbers, etc in 3D data structure. list beh=current, time=, subj=, type=, .......
# TODO sort. Try do.call() stdlib function.
.getAllContexts = function(behaviorName, data, k=NULL, kbefore=0, kafter=0) {
	if (!is.null(k)) {
		kbefore = k;
		kafter = k;
	}
	indicesVector = which(data$behavior == behaviorName);
	indicesVector = indicesVector[indicesVector > kbefore & indicesVector <= (length(data$behavior) - kafter)];
	if (length(indicesVector) == 0) stop(paste("Error: Zero occurances of", behaviorName, "within margins in data provided to .getAllContexts"));
	if (length(indicesVector) < 3) warning(paste("Only", length(indicesVector), "occurances of", behaviorName, "in data provided to .getAllContexts"));
	df = data.frame(centerTimeNum = data$time[indicesVector]);
	for (i in (-kbefore):kafter) {
		df = cbind(df, data$behavior[indicesVector + i]);
		dimnames(df)[[2]][length(dimnames(df)[[2]])] <- paste("n", if(i>=0){"+"}else{""}, as.character(i), sep="");
	}
	print(df);
	centerCol = which(dimnames(df)[[2]] == "n+0");
	# df = df[order(df[,centerCol + 1], df[,centerCol + 2], df[,centerCol + 3]),]      but tailored to actual num of cols.
	return(df);
}









