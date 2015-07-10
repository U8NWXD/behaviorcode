


op = "__June/out_071020141352/"
.calcBasicStats(practice_data, op)
.compareTransitionalProbabilities(practice_data, byTotal = F, paste(op, "byTotalFalse", sep = '_'))
.compareTransitionalProbabilities(practice_data, byTotal = T, paste(op, "byTotalTrue", sep = '_'))
.compareEntropy(practice_data, op)
.makeGroupDotPlots(practice_data, paste(op, "plainMarkov_0.02", sep = '_'), minValForLine = .02)
.makeMulticolorRasterPlots(practice_data, op, colorkey_nostops)
.behavioralDensityGraphs(practice_data, colorkey_startstop, op)