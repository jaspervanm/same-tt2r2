source("scripts/init.R")

study_populations <- list(
	# Information about the distribution of SAMe-TT2R2 scores and TTR (mean+sd)
	# Sometimes the scores aren't completely correct when authors bundle scores (e.g. 0)
	# This is not a problem as long as they don't combine over cutoffs
	Abumuaileq = data.frame(
		# DOI: 10.1093/europace/euu353
		samett2r2_score = 0:4,
		N               = c(247, 425, 174,  46, 19),
		Mean            = c( 62,  58,  54,  55, 53),
		Sd              = c( 18,  18,  19,  15, 18)
	),
	Apostolakis = data.frame(
		# DOI: 10.1378/chest.13-0054
		samett2r2_score = 0:6,
		N               = c( 51, 60, 59, 54, 43, 15,  4),
		Mean            = c( 70, 66, 66, 65, 64, 56, 66),
		Sd              = c( 13, 17, 16, 17, 15, 19, 27)
	),
	Bernaitis = data.frame(
		# DOI: 10.1016/j.jstrokecerebrovasdis.2016.08.030
		samett2r2_score = 2:6,
		N               = c( 339, 517, 210, 56, 15 ),
		Mean            = c( 63.2, 57.3, 54.3, 51.0, 43.4),
		Sd              = c( 34.1, 34.0, 34.7, 32.1, 33.4 )
	),
	Chan = data.frame(
		# DOI: 10.1371/journal.pone.0150674
		# N.B.: very low TTR, minimum score = 2
		samett2r2_score = 2:4,
		N               = c(  254,  646, 442 + 75 + 11),
		Mean            = c( 41.0, 39.0, 36.0 ),
		Sd              = c( 23.3, 25.0, 24.0 )
	),
	Demelo = data.frame(
		# DOI: 10.1016/j.thromres.2016.09.021
		samett2r2_score = c( 0   , 2), # we do not have SDs for separate values
		N               = c( 91  , 44),
		Mean            = c( 64.7, 66.0),
		Sd              = c( 19.5, 20.5)
	),
	Gallego = data.frame(
		# DOI: 10.1016/j.amjmed.2014.05.023
		samett2r2_score = 1:3,
		N               = c( 431   , 332   , 208 ),
		Mean            = c(  79.67,  78.40,  74.25 ),
		Sd              = c(  19.46,  20.28,  20.24 )
	),
	Gorzelak = data.frame(
		# DOI: 10.20452/pamw.3475
		# algebraically derived
		samett2r2_score = 1:2, # i.e. 0-1 and >= 2
		N               = c(7, 97),
		Mean            = c( 64, 50),
		Sd              = c( 26, 28)
	),
	Lip = data.frame(
		# DOI: 10.1016/j.amjcard.2016.08.047
		# N algebraically derived
		samett2r2_score = c( 1, 2), # i.e. 0-1 and >= 2
		N               = c( 119, 110 ),
		Mean            = c( 57.1, 49.8),
		Sd              = c( 22  , 24 )
	),
	Palareti = data.frame(
		# DOI: 10.1160/TH15-10-0830
		samett2r2_score = 0:4,
		N               = c(292, 624, 270, 75, 47),
		Mean            = c( 62,  61,  56, 56, 60),
		Sd              = c( 22,  22,  24, 22, 22)
	),
	Park = data.frame(
		# DOI: 10.5853/jos.2015.17.2.192
		samett2r2_score = 2:6,
		N               = c(  49,  178,  113,   33,  7),
		Mean            = c(58.4, 57.1, 57.8, 55.9, 66.7),
		Sd              = c(17.5, 16.3, 18.6, 17.7, 17.2)
	),
	Poli = data.frame(
		# DOI: 10.1007/s11739-014-1065-8
		samett2r2_score = 0:6,
		N               = c( 189,   435,  288,  127,   44,    5, 1),
		Mean            = c( 71.9, 72.4, 72.0, 68.5, 68.4, 62.6, 52.0 ),
		Sd              = c( 16.9, 14.6, 15.6, 16.9, 15.7, 13.5, 0)
	),
	Roldan = data.frame(
		# DOI: 10.1016/j.amjmed.2015.05.036
		samett2r2_score = 1:2,
		N               = round( 459 * c(0.55, 0.45)),
		Mean            = c(67, 61),
		Sd              = c(18, 16)
	),
	Ruiz = data.frame(
		# DOI: 10.1160/TH15-02-0169
		samett2r2_score = 0:4,
		N               = c( 202,  411,  316,   90, 37),
		Mean            = c(67.5, 64.7, 61.8, 63.0, 52.7),
		Sd              = c(24.6, 26.9, 25.6, 23.2, 28.7)
	)
)

if(!exists("study_findings")) {
	# Read further information about the studies
	# We recalculated test statistics so that
	# SAMe-TT2R2 score >= score_cutoff predicts TTR < qoac_cutoff
	study_findings <- read_excel("search/rapportage_studies_recalculated.xlsx") %>%
		mutate(
			# Create BibTeX reference code for inclusion in manuscript
			Reference = paste0("\\cite{", ref, "}")
		) %>%
		# Filter remarks from Excel file
		select( -opmerkingen ) %>%
		within({
			score_cutoff_included <- ifelse(
				score_cutoff_included == "TRUE" | score_cutoff_included == TRUE | score_cutoff_included == 1,
				TRUE,
				FALSE
			)
			
			# algebraically derive PSEP
			PPV_se <- ( 0.5*(ppv_pe - ppv_ll) + 0.5*(ppv_ul - ppv_pe) ) / 1.96
			NPV_se <- ( 0.5*(npv_pe - npv_ll) + 0.5*(npv_ul - npv_pe) ) / 1.96
			
			PSEP_pe <- ppv_pe + npv_pe - 1
			PSEP_se <- sqrt( PPV_se^2 + NPV_se^2 )
			PSEP_ll <- PSEP_pe - 1.96*PSEP_se
			PSEP_ul <- PSEP_pe + 1.96*PSEP_se
			
			PPV_se <- NPV_se <- PSEP_se <- NULL
			
			# we need to log likelihood ratios
			LRplog_pe <- log(LRp_pe)
			LRplog_ll <- log(LRp_ll)
			LRplog_ul <- log(LRp_ul)
			LRnlog_pe <- log(LRn_pe)
			LRnlog_ll <- log(LRn_ll)
			LRnlog_ul <- log(LRn_ul)
			LRp_pe <- LRp_ll <- LRp_ul <- LRn_pe <- LRn_ll <- LRn_ul <- NULL
		}) %>%
		gather("key", "value", -(study:score_cutoff_included), -Reference ) %>%
		separate("key", c("variable", "attr") ) %>%
		spread( attr, value ) %>%
		within({
			# transform A > X into A >= (X+1)
			# because we use the latter in our computations
			# this is later undone to print the table with study characteristics
			x <- ( score_cutoff_included == FALSE & !is.na( score_cutoff_included) )
			score_cutoff[ x] <- score_cutoff[ x ] + 1
			x <- NULL
			
			# set origin of data
			origin <- "original_data"
			
			# calculate standard error from confidence interval
			se <- ( 0.5*(pe - ll) + 0.5*(ul - pe) ) / 1.96
		})
	
	# Make a supplementary table with study outcomes
	# Again, SAMe-TT2R2 >= cutoff predicts TTR < cutoff
	# This can differ from what the studies originally did, hence other numbers
	# However the original cutoffs are used
	study_findings %>%
		filter( !is.na(score_cutoff), !is.na(qoac_cutoff) ) %>%
		mutate(
			qoac_cutoff  = ifelse(
				is.na( qoac_cutoff ),
				"-",
				paste0("<", qoac_cutoff)
			),
			score_cutoff = ifelse(
				score_cutoff_included == TRUE,
				paste0(">=", score_cutoff),
				paste0(">", score_cutoff - 1) # -1 to undo the +1 from above
			)
		) %>%
		rename( `SAMe-TT2R2` = score_cutoff, TTR = qoac_cutoff ) %>%
		within({
			# undo the log-transformation for LR
			pe[variable == "LRplog"] <- exp(pe[variable == "LRplog"])
			ll[variable == "LRplog"] <- exp(ll[variable == "LRplog"])
			ul[variable == "LRplog"] <- exp(ul[variable == "LRplog"])
			variable[variable == "LRplog"] <- "LRp"
			pe[variable == "LRnlog"] <- exp(pe[variable == "LRnlog"])
			ll[variable == "LRnlog"] <- exp(ll[variable == "LRnlog"])
			ul[variable == "LRnlog"] <- exp(ul[variable == "LRnlog"])
			variable[variable == "LRnlog"] <- "LRn"
		}) %>%
		mutate(
			variable = Variable_names[ variable ],
			pe = round(pe, digits = 2) %>% format( nsmall = 2),
			ll = round(ll, digits = 2) %>% format( nsmall = 2),
			ul = round(ul, digits = 2) %>% format( nsmall = 2),
			report = paste0( pe, " (", ll, "-", ul, ")" ),
			Study = paste0(study, Reference)
		) %>%
		select( Study, `SAMe-TT2R2`, TTR, variable, report ) %>%
		spread( variable, report ) %>%
		do({
			# write the tables to file
			mutate( . ,
					`SAMe-TT2R2` = LaTeX_escape_math(`SAMe-TT2R2`),
					TTR = LaTeX_escape_math(TTR)
			) %>%
				as.data.frame() %>%
				Hmisc::latex(
					title = "Results of individual studies",
					file = "tables/study_findings.tex",
					booktabs = TRUE, table.env = FALSE, rowname = NULL, multicol = FALSE,
					center="centering"
				)

			pandoc.table.return( as.data.frame( . ) ,
								 caption = "Table S2: Results of individual studies  \nLR-, LR+: negative and positive likelihood ratio, respectively; prev: prevalence; PSEP: power of separation; TTR: time in therapeutic range (here also percentage of INR's in therapeutic range)",
								 style = "rmarkdown", split.tables = Inf) %>%
				gsub( ">=", "&ge;", . ) %>%
				gsub( "<=", "&le;", . ) %>%
				gsub( ">", "&gt;", . ) %>%
				gsub( "<", "&lt;", . ) %>%
				write( "tables/study_findings.md" )
			.
		})
	
	study_findings %>%
		# In the Excel file we chose to report statistics in a way that:
		# `SAMe-TT2R2 >= score_cutoff predicts TTR < qoac_cutoff`
		# below we report characteristics the same way the authors did
		within({
			TTR          <- character()
			Score <- character()
			
			hslq <- which( original_study_method == "high score predicts low TTR") # High Score Low Qoac
			hshq <- which( original_study_method == "high score predicts high TTR")
			lshq <- which( original_study_method == "low score predicts high TTR")
			cont <- ( original_study_method %in% c("only continuous TTR", "only in table") )
			
			TTR[ hslq ] <- paste0("<", qoac_cutoff[ hslq ])
			TTR[ c( lshq, hshq) ] <- paste0(">", qoac_cutoff[ c( lshq, hshq) ])
			TTR[ is.na( qoac_cutoff ) | cont] <- "-"
			
			Score[ c( hslq, hshq) ] <- ifelse( score_cutoff_included[ c( hslq, hshq) ] == TRUE,
											paste0(">=", score_cutoff[ c( hslq, hshq) ]),
											paste0(">", score_cutoff[ c( hslq, hshq) ] - 1) ) # to undo the +1 from above
			Score[ lshq ] <- ifelse( score_cutoff_included[ lshq ] == TRUE,
											paste0("<", score_cutoff[ lshq ] ),
											paste0("<=", score_cutoff[ lshq ] -1 ) ) # to undo the +1 from above
			Score[ is.na( score_cutoff ) | cont ] <- "-"
			
			hslq <- hshq <- lshq <- cont <- NULL
		}) %>%
		mutate(
			Study = paste0(study, Reference)
		) %>%
		distinct( Study, Score, TTR, indication, .keep_all = TRUE ) %>%
		select( Study,
				Score, TTR,
				Indication = indication,
				N, Cohort = cohort, `Period excluded` = `period excluded`, `TTR duration`, `TTR method` = ttr_method) %>%
		do({
			# write the tables to file
			mutate( . ,
					Score = LaTeX_escape_math(Score),
					TTR = LaTeX_escape_math(TTR)
				) %>%
				rename( Ind = Indication ) %>%
				as.data.frame() %>%
				Hmisc::latex( title = "Study characteristics",
							  file = "tables/study_characteristics.tex",
							  booktabs = TRUE, table.env = FALSE,
							  rowname = NULL, multicol = FALSE, center="centering")

			pandoc.table.return( as.data.frame( . ) ,
								 caption = "Table 1: Study characteristics  \nAF: Atrial fibrillation; N: number of patients included; period excluded: period excluded in calculation of the TTR; TTR: time in therapeutic range (here also percentage of INR's in therapeutic range); VTE: venous thromboembolism",
								 style = "rmarkdown", split.tables = Inf ) %>%
				gsub( ">=", "&ge;", . ) %>%
				gsub( "<=", "&le;", . ) %>%
				gsub( ">", "&gt;", . ) %>%
				gsub( "<", "&lt;", . ) %>%
				write( "tables/study_characteristics.md")
			.
		})
	
	cutoffs <- study_findings %>%
		# for use in manuscript.R, not in publication
		# do not change the variable names
		distinct( study, score_cutoff, qoac_cutoff, .keep_all = TRUE) %>%
		group_by( score_cutoff, qoac_cutoff) %>%
		summarise(
			Npat = sum(N),
			Nstud = n(),
			refs = paste0("\\cite{", paste(ref, sep = "", collapse = ", " ), "}")
		) %>%
		arrange( score_cutoff, qoac_cutoff ) %>%
		ungroup() %>%
		do({
			# write the tables to file
			pandoc.table.return( as.data.frame( . ) , style = "rmarkdown", split.tables = Inf) %>% write( "tables/cutoffs.md" )
			.
		})
}
