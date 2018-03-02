# MANUSCRIPT.R
# PURPOSE: SET VARIABLES FOR USE IN THE MANUSCRIPT

source("scripts/compare_studies_simulations.R")
# : calls init
# : calls simulation: creates simulation / simulation_summary / simulation_and_studies
# : calls studies: writes tables/studies
# : writes plots/correlation_simulation_original.png
# : creates studies_simulation_stats (correlation, mean and sd from difference)
# : creates outside_CI_mean

source("scripts/meta_analysis.R")
# : creates MetaAnalysis
# : writes tables/meta_analysis
# : writes plots/forest_plots

source("scripts/plot_pretest_posttest.R")
# : writes plots/pretest_posttest

source("scripts/score_effect_on_TTR.R")
# : writes plot/forest_plot_ttr
# : creates ttr_comp with standardised mean difference in TTR

source("scripts/sensitivity_analyses.R")
# : writes plots/sensitivity_analysis_*
# : writes tables/sensitivity_analysis

# General Functions
Range <- function( Variable, grouping_variable = NULL, Data = study_findings) {
	if( !is.null( grouping_variable ) ) {
		Data <- Data %>%
			group_by_( as.name( grouping_variable) )
	}
	
	Data %>%
		filter( variable == Variable ) %>%
		summarise(
			min = min( pe, na.rm = TRUE ),
			max = max( pe, na.rm = TRUE)
		) %>%
		mutate(
			range = paste(
				round( min, digits = 2),
				round( max, digits = 2),
				sep = "-" ),
			pct_range = map2_chr( min, max, ~Pct( c( .x, .y ) ) )
		)
}
paste_refs <- function(x) {
	if(is.data.frame(x)) {
		x <- x %>%
			distinct(ref) %>%
			unlist()
	}
	x %>%
		paste(collapse=", ") %>%
		paste0("\\cite{", . , "}")
}

# Introduction


# Methods


# Results: articles identified
N <- list(
	studies = study_findings %>% distinct(study) %>% unlist() %>% length(),
	indication = study_findings %>% distinct(study, indication) %>% with( table( indication) ) %>% c()
)
ind_refs <- study_findings %>%
	distinct(ref, .keep_all = TRUE) %>%
	with( tapply( ref, indication, paste_refs ))
coh_refs <- study_findings %>%
	distinct(ref, .keep_all = TRUE) %>%
	with( tapply( ref, cohort, paste_refs ))
coh_nr_refs <- study_findings %>%
	filter( cohort %in% c("not reported", "mixed") ) %>%
	paste_refs()

Nmeta <- simulation_and_studies %>%
	filter(	score_cutoff %in% 2:3,
			qoac_cutoff == 70,
			se > 0,
			variable == "sens" # dit haalt extra statistieken uit de simulatie eruit
	) %>%
	distinct( study ) %>%
	unlist() %>%
	length()
origin <- simulation_and_studies %>%
	filter(	score_cutoff %in% 2:3,
			qoac_cutoff == 70,
			se > 0,
			variable == "LRplog" # dit haalt extra statistieken uit de simulatie eruit
	) %>%
	distinct( study, score_cutoff, qoac_cutoff, .keep_all = TRUE ) %>%
	group_by( qoac_cutoff, score_cutoff ) %>%
	do(
		table(.$origin) %>% as.data.frame()
	) %>%
	spread( Var1, Freq) %>%
	mutate(
		total = original_data + simulation
	)

# Results: individual studies
most_used_cutoff <- cutoffs %>%
	arrange( desc(Nstud ) ) %>%
	slice( 1 ) %>%
	unlist()
no_cutoff <- cutoffs %>%
	filter( is.na(score_cutoff), is.na(qoac_cutoff)) %>%
	slice( 1 ) %>%
	unlist()
high_score_range <- study_findings %>%
	filter( !is.na( score_cutoff ), !is.na( qoac_cutoff ) ) %>%
	group_by( study, score_cutoff, qoac_cutoff) %>%
	select( variable, pe, study, score_cutoff, qoac_cutoff ) %>%
	spread( variable, pe) %>%
	summarise(
		score_above_cutoff = sens * prev + (1 - spec) * (1 - prev)
	) %>%
	group_by( score_cutoff ) %>%
	summarise(
		p_high_score = Pct( range( score_above_cutoff ) )
	) %>%
	filter( score_cutoff %in% 2:3 )
cohort <- study_findings %>%
	distinct( study, .keep_all = TRUE ) %>%
	with( table( cohort ) )
unexpected_LR <- study_findings %>%
	filter(
		(variable == "LRplog" & exp(pe) < 1) |
			(variable == "LRnlog" & exp(pe) > 1)
	) %>%
	distinct( study, ref ) %>%
	mutate( ref = as.character(ref))

# Results: validation of the simulation
SimVal <- c(
	outside_CI = Pct( 1 - outside_CI_mean ),
	Cor  = Pct( studies_simulation_stats["cor"]),
	ceiling( studies_simulation_stats["mean"] * 100) / 100, # in abstract we say "less than X"
	round( studies_simulation_stats["sd"], digits = 2)
)

# Results: meta-analysis
I2_sens_spec <- MetaAnalysis %>%
	filter(variable %in% c("sens", "spec")) %>%
	select(score_cutoff, variable, qoac_cutoff, I2) %>%
	unnest() %>%
	with(min(I2)) %>%
	trunc()

MRange <- function( Variable, Score_cutoff = 2:4, Data = MetaAnalysis ) {
	Data %>%
		filter( variable == Variable, score_cutoff %in% Score_cutoff ) %>%
		mutate(
			range = paste(
				format( round( ll, digits = 2), nsmall = 2),
				format( round( ul, digits = 2), nsmall = 2),
				sep = "-" ),
			pct_range = map2_chr( ll, ul, ~Pct( c( .x, .y ) ) )
		)
}
MPct <- function( Data, CI = TRUE, CItext = FALSE ) {
	if( nrow( Data ) > 1) {
		alply( Data, 1, MPct, CI )
	} else {
		paste0(
			Pct( Data$pe ),
			ifelse( CI == FALSE,
					"",
					paste0(
						ifelse( CItext == TRUE,
								" (95% CI ",
								" (" ),
						Data$pct_range, ")")
			)
		)
	}
}
pe_ci <- function( ... ) {
	MRange( ... ) %>%
		alply( 1, function( x ) {
			paste0( round(x$pe, 2),
					" (", x$range, ")")
		}) %>%
		unlist()
}

# Discussion
cutoff_use <- MetaAnalysis %>%
	filter( variable %in% c("LRn", "LRp")) %>%
	select( score_cutoff, variable, pe ) %>%
	spread( variable, pe) %>%
	mutate(
		LRp_pop = do_LR( 1 / LRp , 0.70 ) %>% round( 3 ) * 100,
		LRn_pop = do_LR( 1 / LRn , 0.20 ) %>% round( 3 ) * 100
	)
