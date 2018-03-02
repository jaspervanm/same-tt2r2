source("scripts/init.R")

# start with function definitions,
# the actual simulation is started further below

# file created to save the results of the computation
# this file is loaded if it exists, and then the computation will be skipped
generated_datafile <- "data/simulated_study_performance.Rdata"

# estBetaParams function via StackOverflow, user "Max"
# URL http://stats.stackexchange.com/questions/12232
# function to construct a beta distribution based on parameters from normal distribution
# TTR cannot be a normal distribution, because it ranges from 0 to 1
estBetaParams <- function(mu, var) {
	alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
	beta  <- alpha * (1 / mu - 1)
	return(params = list(
		alpha = alpha,
		beta = beta) )
}

generate_sample <- function(population_data) {
	# here we use the study data to simulate a study population
	# For example:
	# create 100 people with score X and TTR mean Y, sd Z
	# because TTR cannot be normally distributed (has to be in range 0-100%)
	# we use a beta distribution (parameter estimation above)
	adply( population_data, 1, function(x) {
		betaP <- with(x, estBetaParams( Mean/100, (Sd/100)^2) )
		data.frame(
			samett2r2_score = x$samett2r2_score,
			value           = rbeta( x$N, betaP$alpha, betaP$beta)*100
		)
	} )
}

simulate_samett2r2 <- function(population_data, qoac_cutoff, score_cutoff, iter = 1000 ) {
	# so standard 1000 iterations per study
	
	ldply( seq_len( iter ), function( i ) {
		#{ building the population
		population <- generate_sample( population_data )
		#} building the population
		
		#{ we will assess the following cutoffs
		cutoffs <- expand.grid(
			qoac_cutoff  = qoac_cutoff,
			score_cutoff = score_cutoff
		)
		
		adply( cutoffs, 1, function( cutoff ) {
			
			#{ creating intermediary objects
			population$reference = ( population$value < cutoff$qoac_cutoff )
			population$predictor = ( population$samett2r2_score >= cutoff$score_cutoff )
			
			test.tab <- population %>%
				mutate(
					pred = factor( predictor, levels = c( TRUE, FALSE) ),
					ref  = factor( reference, levels = c( TRUE, FALSE) )
				) %>%
				with( table(pred, ref) )
			
			prev.ref  <- mean( population$reference )
			prev.pred <- mean( population$predictor )
			#} creating intermediary objects
			
			#{ returning final object
			output <- data.frame(
				indication   = ifelse( !is.null( population_data[1, "indication"] ),
									   as.character( population_data[1, "indication"] ),
									   NA ),
				N            = nrow(population),
				qoac_cutoff  = cutoff$qoac,
				score_cutoff = cutoff$score,
				prev.ref     = prev.ref,
				prev.pred    = prev.pred,
				sens         = sensitivity( test.tab ),
				spec         = specificity( test.tab ),
				npv          = negPredValue( test.tab ),
				ppv          = posPredValue( test.tab ),
				PSEP         = posPredValue( test.tab ) + negPredValue( test.tab ) - 1, # Altman, PMID 10694730
				LRplog       = log( sensitivity( test.tab ) / (1 - specificity( test.tab ) ) ),
				LRnlog       = log( (1 - sensitivity( test.tab ) ) / specificity( test.tab ) )
			)
			output[ 1, sapply( output, is.infinite) ] <- NaN
			output
			#} returning final object
		})
		
	})
}

# load information about the studies
source("scripts/studies.R")

# set the parameters
# the parameters are broader than our research question
# to allow us to check how well the simulation performs
params <- lapply( study_populations, list,
				  qoac_cutoff  = c(65, 70),
				  score_cutoff = c(2, 3, 4) )
# make study names an explicit variable for further reference
for(i in seq_along(params)) {
	params[i][[1]][[1]]$study <- names(params)[i]
}

# The actual simulation starts below
# But first check if it actually needs to be run, because it takes long
if( exists( "simulation") ) {
	# nothing to do,
	# script was probably called again by other script.
} else if( file.exists(generated_datafile) ) {
	# load cache if it exists
	load(generated_datafile)
} else  {
	# we fix the random number generator,
	# otherwise the results will differ slightly every time the simulation is run
	set.seed(123)
	
	# do the actual simulation
	simulation <- ldply( params, do.call, what = simulate_samett2r2,
						 .id = "study" )
	
	# save results to save time next time
	save( simulation, file = generated_datafile )
}

# Summarise the 1000 simulations for one study
simulation_summary <- simulation %>%
	mutate(
		study = as.character(study)
	) %>%
	group_by(study, qoac_cutoff, score_cutoff, indication, N) %>%
	summarise_all( c("mean", "var", "RIlow", "RIhigh"), na.rm = TRUE ) %>%
	gather( "key", "value", -study, -qoac_cutoff, -score_cutoff, -indication, -N ) %>%
	separate( key, c("variable", "stat"), "_") %>%
	within({
		stat[stat == "mean"]   <- "pe" # point estimate
		stat[stat == "RIlow"]  <- "ll" # lower limit
		stat[stat == "RIhigh"] <- "ul" # upper limit
	}) %>%
	spread(stat, value) %>%
	mutate(
		se     = sqrt( var ),
		var    = NULL,
		origin = "simulation"
	)
simulation_summary$variable[ simulation_summary$variable == "prev.ref" ] <- "prev"

# some studies can only be used for certain cutoffs,
# but we simulate them all
# below we remove the simulated cutoffs that we shouldn't use
simulation_summary <- simulation_summary %>%
	filter(
		# leave out the ones where no patients with score < 2
		!( study == "Bernaitis" & score_cutoff == 2),
		!( study == "Chan" & score_cutoff == 2),
		!( study == "Park" & score_cutoff == 2),
		# leave out the ones with insufficient data to evaluate score >= 3
		!( study == "Demelo" & score_cutoff == 3),
		!( study == "Gorzelak" & score_cutoff == 3),
		!( study == "Lip" & score_cutoff == 3),
		!( study == "Lobos" & score_cutoff == 3),
		!( study == "Roldan" & score_cutoff == 3),
		!( study == "Szymanski" & score_cutoff == 3)
	)

# we determine mean and standard deviation of one sample of generated TTR's
# for sensitivity analyses
set.seed(123)
qoac_stats <- map_df(names(params), function(j) {
	x <- generate_sample(params[[j]][[1]])
	
	cbind(
		simulation_summary %>%
			ungroup() %>%
			filter( study == j, variable == "prev") %>%
			distinct(study, qoac_cutoff, score_cutoff, .keep_all = TRUE) %>%
			select( study, qoac_cutoff, score_cutoff, prev_low_ttr = pe),
		mean_ttr = mean(x$value),
		sd_ttr   = sd(x$value)
	)
}) %>%
	group_by(score_cutoff, qoac_cutoff) %>%
	mutate(
		mean_ttr_cat = cut(mean_ttr, c(0, median(mean_ttr, na.rm = TRUE), 100), labels = c("TTR low", "TTR high"))	,
		sd_ttr_cat = cut(sd_ttr, c(0, median(sd_ttr, na.rm = TRUE), Inf), labels = c("SD low", "SD high")),
		qoac_cat = paste(mean_ttr_cat, sd_ttr_cat, sep = ", "),
		prev_low_ttr_cat = cut(prev_low_ttr, c(0, median(prev_low_ttr, na.rm = TRUE), 1), labels = c("few low TTR", "many low TTR") )
	) %>%
	ungroup()

simulation_and_studies <- bind_rows(
	# data from the original studies supersedes simulated data
	semi_join( study_findings, simulation_summary, by = c("study", "qoac_cutoff", "score_cutoff", "variable") ),
	# take all data from the simulation where no original data exists
	anti_join( simulation_summary, study_findings, by = c("study", "qoac_cutoff", "score_cutoff", "variable") ),
	# take all data from the original study where no simulation exists
	# e.g. when no parameters for normal distribution were provided
	anti_join( study_findings, simulation_summary, by = c("study", "qoac_cutoff", "score_cutoff", "variable") )
) %>%
	select(-year, -ref, -indication, -ttr_method) %>% # deel is leeg
	left_join( study_findings %>% distinct(study, year, ref, indication, ttr_method), by = "study") %>%
	mutate_at(vars(indication, ttr_method), as.factor) %>%
	left_join( qoac_stats, by = c("study", "qoac_cutoff", "score_cutoff"))
