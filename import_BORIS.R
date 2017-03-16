# behavior_log = import_BORIS('my_log.csv', 15)

import_BORIS = function(file, start_row_num) {
	input_log = read.csv('~/Downloads/csv test.txt.csv', header = T, skip = start_row_num)
	log = data.frame(
		time = input_log$Time,
		behavior = input_log$Behavior,
		subject = input_log$Subject,
		type = input_log$X,
		pair_time = NA,
		duration = NA
	)
	log$type[log$type == 'START'] <- 'start'
	log$type[log$type == 'STOP'] <- 'stop'
	log$type[log$type == 'POINT'] <- 'neither'
	
	log = .BORISpairStartStop(log)
	
	attr(log, 'assay.start') <- list(mark = NA, time = 0);
	return(log)
}

.BORISpairStartStop = function(log) {
	ssbehs = names(table(log$behavior[log$type %in% c('start', 'stop')]))
	for (beh in ssbehs) {
		start_indices = which(log$behavior == beh & log$type == 'start')
		stop_indices = which(log$behavior == beh & log$type == 'stop')
		if (stop_indices[1] < start_indices[1]) {
			stop_indices = stop_indices[-1]
		}
		if (stop_indices[length(stop_indices)] < start_indices[length(start_indices)]) {
			start_indices = start_indices[-length(start_indices)]
		}
		if (length(start_indices) != length(stop_indices)) {
			stop('Unequal number of starts and stops')
		}
		
		for (i in 1:length(start_indices)) {
			if (start_indices[i] >= stop_indices[i]) {
				stop('Mismatched start and stop')
			}
			log$pair_time[start_indices[i]] <- log$time[stop_indices[i]]
			log$pair_time[stop_indices[i]] <- log$time[start_indices[i]]
			log$duration[c(start_indices[i], stop_indices[i])] = log$time[stop_indices[i]] - log$time[start_indices[i]]
		}
	}
	return(log)
}
