source("scripts/meta_analysis.R")

sensitivity_variables <- c(
	indication       = "Indication",
	ttr_method       = "TTR method",
	origin           = "Origin",
	prev_low_ttr_cat = "Prevalence of low TTR",
	mean_ttr_cat     = "Categorised mean TTR",
	sd_ttr_cat       = "Categorised sd TTR",
	qoac_cat         = "Quality of anticoagulation"
)

Sensitivity_Analysis <- map_df( names(sensitivity_variables), function(s) {
	MetaAnalysis %>%
		within({
			# all (likelihood) ratios are log-transformed
			lr <- variable %in% c("LRp", "LRn")
		}) %>%
		group_by(score_cutoff, qoac_cutoff, variable, lr) %>%
		with({
			pmap_df(
				list(data, score_cutoff, qoac_cutoff, variable),
				function(Data, score_cutoff, qoac_cutoff, variable) {
					Data %>%
						group_by_( s ) %>%
						do({
							with( . , {
								x <- rma(pe, sei = se, weights = N,
										 measure = ifelse( all( lr ), "RR", "GEN" ),
										 control = list( stepadj = 0.25) )

								data.frame(
									sensitivity_variable = sensitivity_variables[ s ],
									sensitivity_value = as.character( unlist(.[1, s]) ),
									score_cutoff = score_cutoff,
									qoac_cutoff = qoac_cutoff,
									variable = variable,
									pe = x$b,
									ll = x$ci.lb,
									ul = x$ci.ub,
									p  = x$pval,
									N  = nrow( . ),
									stringsAsFactors = FALSE
								)
							})
						})
				})
		}) %>%
		within({
			sensitivity_value[sensitivity_value == "original_data"] <- "original data"
			lr <- variable %in% c("LRp", "LRn")
			pe[ lr ] <- exp( pe[ lr ])
			ll[ lr ] <- exp( ll[ lr ])
			ul[ lr ] <- exp( ul[ lr ])
			plot_legend <- paste0( sensitivity_variable, ": ", sensitivity_value)
			rm(lr)
		})
}) %>%
	rbind.fill(
		MetaAnalysis %>%
			select(-data, -meta_analysis) %>%
			mutate(
				sensitivity_variable = "Full analysis",
				plot_legend = "Full analysis",
				sensitivity_value = ""
			)
	) %>%
	mutate(
		sensitivity_variable = factor(sensitivity_variable, levels = c("Full analysis", sensitivity_variables))
	)

# Summary table
Sensitivity_Analysis %>%
	filter(
		variable %in% c("LRp", "LRn", "PSEP"),
		!is.na(sensitivity_value)
	) %>%
	mutate(
		pe       = format( pe, digits = 0, nsmall = 2),
		ll       = format( ll, digits = 0, nsmall = 2),
		ul       = format( ul, digits = 0, nsmall = 2),
		ci       = paste0( "(", ll, "-", ul, ")" ),
		value    = paste(  pe, ci),
		cutoff   = paste(">=", score_cutoff),
		variable = Variable_names[ variable ],
		`Sensitivity analysis` = plot_legend,
		N        = as.integer( N )
	) %>%
	group_by(`Sensitivity analysis`, cutoff, N) %>%
	select( variable, value, contains("sensitivity") ) %>%
	spread( variable, value ) %>%
	arrange( cutoff, sensitivity_variable, sensitivity_value) %>%
	ungroup() %>%
	select(cutoff, sensitivity_variable, sensitivity_value, N, contains("LR"), PSEP) %>%
	group_by(cutoff) %>%
	do({
		x <- .
		pandoc.table.return( as.data.frame( x ) ,
							 caption = "S7 Table: Sensitivity analyses.  \n
LR-, LR+: negative and positive likelihood ratio, respectively; PSEP: power of separation.\n
AF: atrial fibrillation; N: number of studies; PINRR: proportion of INR’s in range; TTR: time within therapeutic range; VTE: venous thrombo-embolism.\n
*\tCategorisation of prevalence of low TTR, mean TTR and standard deviation TTR based on upper and lower half.\n
**\tQuality of anticoagulation is a combination of categories from mean TTR and sd TTR.
",
	 						 style = "rmarkdown", split.tables = Inf ) %>%
	 		write( paste0("tables/sensitivity_analysis_", substr(x[["cutoff"]][1], 4, 4), ".md" ) )
	 	x
	 })


do_sensitivity_plot <- function(Score_cutoff) {
	Sensitivity_Analysis %>%
		filter(
			variable %in% c("sens", "spec", "PSEP"),
			score_cutoff == Score_cutoff,
			!is.na(sensitivity_value)
		) %>%
		ggplot( aes(y = pe, ymin = ll, ymax = ul, colour = sensitivity_variable,
					x = paste(sensitivity_value, paste0("(", N, ifelse(N == 1, " study)", " studies)") ) ) ) ) +
		geom_point( position = position_dodge(width = 0.5)) +
		geom_linerange(	position = position_dodge(width = 0.5) ) +
		facet_grid( sensitivity_variable ~ Variable_names[variable], scales = "free", space = "free", switch = "y" ) +
		theme_few() +
		theme(
			axis.text.y = element_text(hjust = 0),
			legend.position = "bottom",
			strip.placement = "outside",
			strip.text.y = element_text(angle = 180, vjust = 0.85, hjust = 0)
		) +
		scale_colour_manual(values = unname(colours), name = NULL, guide = FALSE) +
		scale_linetype_identity() +
		labs(x = NULL, y = "estimate (95% CI)", title = paste0("Performance of SAMe-TT2R2 >= ", Score_cutoff) ) +
		coord_flip(ylim = 0:1)
}
sensitivity_plot <- arrangeGrob(
	do_sensitivity_plot(2),
	do_sensitivity_plot(3),
	ncol = 1
)
ggsave( "sensitivity_analysis_02.png", sensitivity_plot, path = "plots/", width = 1000/DPI, height = 1000/DPI, dpi = DPI)
ggsave( "sensitivity_analysis_02.pdf", sensitivity_plot,path = "plots/", width = 1000/DPI, height = 1000/DPI, dpi = DPI)


do_LR_sensitivity_plot <- function(Variable, Score_cutoff, breaks, limits = range(breaks)) {
	p <- ggplot(
		Sensitivity_Analysis %>%
			filter(
				score_cutoff == Score_cutoff,
				variable == Variable,
				!is.na(sensitivity_value)
			),
		aes( x = paste(sensitivity_value, paste0("(", N, ifelse(N == 1, " study)", " studies)") ) ),
			 y = pe, ymin = ll, ymax = ul,
			 colour = sensitivity_variable,
			 linetype = ifelse(ll < limits[1] | ul > limits[2], "dashed", "solid")) ) +
		# Add reference line
		geom_hline( aes( yintercept = 1 ), linetype = "dotted", colour = "gray" ) +
		geom_linerange( ) +
		geom_point() +
		facet_grid(
			sensitivity_variable ~ paste(Variable_names[Variable], "of SAMe-TT2R2 >=", score_cutoff),
			scales = "free", space = "free", switch = "y"
		) +
		theme_few() +
		theme(
			strip.placement = "outside",
			strip.text.y = element_text(angle = 180, vjust = 0.85, hjust = 0)
		) +
		scale_x_discrete( ) +
		scale_y_log10(breaks = breaks) +
		scale_colour_manual(values = unname(colours), guide = FALSE) +
		scale_linetype_identity(guide = FALSE) +
		xlab("") + ylab("estimate (95% CI)") +
		coord_flip(ylim = limits)
	
	if (Variable == "LRp") {
		p <- p + theme(axis.title.y = element_blank(),
					   axis.text.y  = element_blank(),
					   axis.ticks.y = element_blank(),
					   strip.background = element_blank(),
					   strip.text.y = element_blank()
		)
	}
	p
}

LR_sensitivity_plot <- arrangeGrob(
	do_LR_sensitivity_plot("LRn", 2, breaks = c(0.67, 1, 1.42)),
	do_LR_sensitivity_plot("LRp", 2, breaks = c(0.8, 1, 1.5, 2.25, 3)),
	do_LR_sensitivity_plot("LRn", 3, breaks = c(0.67, 1, 1.42)),
	do_LR_sensitivity_plot("LRp", 3, breaks = c(0.8, 1, 1.5, 2.25, 3)),
	ncol = 2,
	widths = c(4.5,2)
)

ggsave( "sensitivity_analysis_01.png", LR_sensitivity_plot, path = "plots/", width = 950/DPI, height = 1050/DPI, dpi = DPI)
ggsave( "sensitivity_analysis_01.pdf", LR_sensitivity_plot, path = "plots/", width = 9.5, height = 10.5, dpi = 450)
