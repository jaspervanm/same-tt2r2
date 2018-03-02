source("scripts/simulation.R")

# CALCULATING THE META-ANALYSIS
MetaAnalysis <- simulation_and_studies %>%
	within({
		# all (likelihood) ratios have been log-transformed
		lr <- variable %in% c("LRplog", "LRnlog")
		slab <- paste0(study, " (", year, ")")
	}) %>%
	# selecting suitable data
	filter(
		score_cutoff %in% c(2, 3), # cutoff for our research question
		qoac_cutoff  == 70,        # cutoff for our research question
		variable != "prev.pred"    # this was given, so se = 0
	) %>%
	# alphabetically list studies for forest plot
	arrange( study ) %>%
	# performing the meta-analysis
	group_by( score_cutoff, variable, qoac_cutoff ) %>%
	nest() %>%
	mutate(
		meta_analysis = map( data, ~with( . , {
			rma(pe, sei = se, weights = N, slab = slab,
				measure = ifelse( all( lr ), "RR", "GEN" ),
				control = list( stepadj = 0.25) )
		} ) ) ,
		# extracting point estimate (pe),
		# lower (ll) and upper limit (ul) of the confidence interval
		pe = map_dbl( meta_analysis, "b"),
		ll = map_dbl( meta_analysis, "ci.lb" ),
		ul = map_dbl( meta_analysis, "ci.ub" ),
		p  = map_dbl( meta_analysis, "pval"),
		N  = map_dbl( data, nrow),
		I2 = map( meta_analysis, ~confint( . )$random[3, ]),
		H2 = map( meta_analysis, ~confint( . )$random[4, ]),
		Q  = map_dbl( meta_analysis, "QE"),
		Qp = map_dbl( meta_analysis, "QEp")
	) %>%
	# undoing the log-transformation for ratios
	within({
		lr <- variable %in% c("LRplog", "LRnlog")
		variable[variable == "LRplog"] <- "LRp"
		variable[variable == "LRnlog"] <- "LRn"
		pe[ lr ] <- exp( pe[ lr ])
		ll[ lr ] <- exp( ll[ lr ])
		ul[ lr ] <- exp( ul[ lr ])
		rm( lr )
	})

# CREATING AND PRINTING THE FOREST PLOTS
ticks <- list(
	LRp =  c( 0.5, 1, 2, 3),
	LRn =  c( 0.5, 1, 1.5),
	sens = seq( 0, 1, by = 0.25),
	spec = seq( 0, 1, by = 0.25),
	PSEP = c( -0.25, 0, 0.25)
)

# General function to create the forest plots
do_forest_plot <- function(Variable, Score_cutoff, breaks, limits = range(breaks), suppressText = FALSE, pool = TRUE) {
	x <- MetaAnalysis %>%
		filter(variable == Variable, score_cutoff == Score_cutoff) %>%
		rename(
			ll_total = ll,
			ul_total = ul,
			pe_total = pe,
			N_total  = N
		)
	
	Data <- x %>%
		unnest(data) %>%
		mutate(
			study_legend = paste0(study, " (", year, ")" )
		) %>%
		left_join(
			influence(x[[1, "meta_analysis"]])$inf %>%
				rownames_to_column("slab") %>%
				select(slab, weight),
			by = "slab"
		)
	I2         <- confint(x[[1, "meta_analysis"]])$random[3, ]
	I2_caption <- paste0(
		"I2 = ",
		round(I2["estimate"], 1),
		"% (95% CI ", round( I2["ci.lb"], 1 ), " - ", round(I2["ci.ub"], 1),
		", p ", ifelse(x[["Qp"]] < 0.001, " < 0.001", paste("=", round(x[["Qp"]], 3) ) ), ")"
	)
	
	if(Variable %in% c("LRp", "LRn")) {
		Data$ll <- exp(Data$ll)
		Data$pe <- exp(Data$pe)
		Data$ul <- exp(Data$ul)
	}
	
	if(pool == TRUE) {
		pooled_data <- data.frame(
			study_legend = "pooled data",
			variable = "pooled",
			score_cutoff = Score_cutoff,
			ll = Data[[1, "ll_total"]],
			ul = Data[[1, "ul_total"]],
			pe = Data[[1, "pe_total"]],
			N  = sum(Data$N),
			stringsAsFactors = FALSE,
			weight = 100
		)
		Data <- plyr::rbind.fill( Data, pooled_data)
		Limits <- c( pooled_data$study_legend, sort(Data[-nrow(Data), "study_legend"], TRUE) )
		
		Caption <- paste0("p (random effects) ",
						  ifelse(x[["p"]] < 0.001, "< 0.001", paste("=", round(x[["p"]], 3))),
						  "\n", I2_caption)
	} else {
		Limits  <- sort(Data$study_legend, TRUE)
		Caption <- I2_caption
	}
	
	p <- ggplot(
		Data,
		aes( x = study_legend,
			 y = pe, ymin = ll, ymax = ul,
			 colour   = variable,
			 linetype = ifelse(ll < limits[1] | ul > limits[2], "dashed", "solid"),
			 shape    = ifelse(!is.na(study), "study", "pooled")
		) ) +
		geom_linerange( ) +
		geom_point( aes( size = round(25 * sqrt( weight / sum(weight)) ) ) ) +
		theme_few(base_size = 10) +
		theme(
			strip.placement = "outside",
			strip.text.y = element_text(angle = 180, vjust = 0.85, hjust = 0)
		) +
		scale_x_discrete( limits = Limits ) +
		scale_colour_manual(values = colours[unique(Data$variable)], guide = FALSE) +
		scale_linetype_identity(guide = FALSE) +
		scale_shape_manual( guide = FALSE, values = c(study = 1, pooled = 5) ) +
		scale_size(guide = FALSE) +
		labs(
			x       = "",
			y       = "estimate (95% CI)",
			title   = paste(Variable_names[Variable], "of SAMe-TT2R2 >=", Score_cutoff),
			caption = Caption
		) +
		coord_flip(ylim = limits) +
		theme(plot.caption = element_text(hjust = 0, size = 8, face = "italic"))
	
	if (Variable %in% c("LRp", "LRn")) {
		p <- p + scale_y_log10(breaks = breaks) +
			# Add reference line
			geom_hline( aes( yintercept = 1 ), linetype = "dotted", colour = "gray" )
	} else {
		p <- p + scale_y_continuous(breaks = breaks)
		if(Variable == "PSEP") {
			p <- p + 
				# Add reference line
				geom_hline( aes( yintercept = 0 ), linetype = "dotted", colour = "gray" )
		}
	}
	
	if (suppressText == TRUE) {
		p <- p + theme(axis.title.y = element_blank(),
					   axis.text.y  = element_blank(),
					   axis.ticks.y = element_blank(),
					   strip.background = element_blank(),
					   strip.text.y = element_blank()
		)
	}
	p
}

# Main forest plot: likelihood ratios
LR_plot <- arrangeGrob(
	do_forest_plot("LRn", 2, breaks = ticks$LRn),
	do_forest_plot("LRp", 2, breaks = ticks$LRp, suppressText = TRUE),
	do_forest_plot("LRn", 3, breaks = ticks$LRn),
	do_forest_plot("LRp", 3, breaks = ticks$LRp, suppressText = TRUE),
	ncol = 2,
	widths = c(3,2)
)
ggsave("plots/forest_plots.png", LR_plot, width = 5.2, height = 5.2, scale = 1.4)
ggsave("plots/forest_plots.eps", LR_plot, width = 5.2, height = 5.2, scale = 1.4)

# Supplementary forest plots: sens, spec, PSEP
sup_plot <- arrangeGrob(
	do_forest_plot("PSEP", 2, breaks = ticks$PSEP),
	do_forest_plot("PSEP", 3, breaks = ticks$PSEP, suppressText = TRUE),
	do_forest_plot("sens", 2, breaks = ticks$sens, pool = FALSE),
	do_forest_plot("sens", 3, breaks = ticks$sens, suppressText = TRUE, pool = FALSE),
	do_forest_plot("spec", 2, breaks = ticks$spec, pool = FALSE),
	do_forest_plot("spec", 3, breaks = ticks$spec, suppressText = TRUE, pool = FALSE),
	ncol = 2,
	widths = c(3, 2)
)
ggsave("plots/forest_plots_sup.png", sup_plot, width = 6, height = 9, scale = 1.3)
ggsave("plots/forest_plots_sup.pdf", sup_plot, width = 6, height = 9, scale = 1.3)

# CREATING AND PRINTING THE FUNNEL PLOTS
LRFunnel <- study_findings %>%
	filter( variable %in% c("LRplog", "LRnlog"),
			is.finite(pe),
			pe != 0 ) %>%
	arrange( study ) %>%
	group_by( score_cutoff, variable, qoac_cutoff ) %>%
	mutate( Nstud = n() ) %>%
	group_by( Nstud, add = TRUE ) %>%
	by_slice(
		with, rma( pe, sei = se, weights = N, slab = study,
				   measure = "RR",
				   control = list( stepadj = 0.25) ), .to = "meta_analysis"
	)

postscript("plots/funnel_plots.eps", onefile = FALSE,
		   width = 11, height = 11,
		   paper = "special", horizontal = FALSE)
par(mfcol = c(2,2)) # combine figures into one

LRFunnel %>%
	filter( Nstud >= 3 ) %>%
	with({
		map2(meta_analysis,
			 paste(
			 	"Funnel plot of log", Variable_names[gsub("log", "", variable)],
			 	"for SAMe-TT2R2 >=", score_cutoff,
			 	"to predict TTR <", qoac_cutoff),
			 ~funnel( .x, xlab = .y ) )
	} ) %>%
	invisible()

par(mfcol = c(1,1))
dev.off()

# TEST FOR FUNNEL PLOT ASYMMETRY
pub_bias <- LRFunnel %>%
	filter( Nstud > 2 ) %>%
	mutate(
		memrm = map( meta_analysis, ~regtest( . ) )
	)
pub_bias$memrm
map_dbl(pub_bias$memrm, c("pval"))

# CREATING A SUMMARY TABLE
MetaAnalysis %>%
	filter( variable %in% c("LRp", "LRn", "PSEP")) %>%
	mutate(
		pe       = format( pe, digits = 0, nsmall = 2),
		ll       = format( ll, digits = 0, nsmall = 2),
		ul       = format( ul, digits = 0, nsmall = 2),
		ci       = paste0( "(", ll, "-", ul, ")" ),
		value    = paste(  pe, ci),
		cutoff   = paste(">=", score_cutoff),
		variable = Variable_names[ variable ]
	) %>%
	select( cutoff, variable, value ) %>%
	spread( variable, value ) %>%
	do({
			mutate( . ,
				cutoff = LaTeX_escape_math(cutoff)
			) %>%
			rename( `SAMe-TT\\textsubscript{2}R\\textsubscript{2}` = cutoff ) %>%
			as.data.frame() %>%
			Hmisc::latex(
				title = "Meta-analysis",
				file = "tables/meta_analysis.tex",
				booktabs = TRUE, table.env = FALSE, rowname = NULL, multicol = FALSE, center="centering"
			)
		
		pandoc.table.return( as.data.frame( . ) ,
							 caption = "Table 2: Performance of the SAMe-TT2R2 score to predict TTR \\\\<70\\%.  \nLR-, LR+: negative and positive likelihood ratio, respectively; PSEP: power of separation; TTR: time in therapeutic range",
							 style = "rmarkdown", split.tables = Inf ) %>%
			write( "tables/meta_analysis.md" )
		.
	})
