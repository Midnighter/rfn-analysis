
library(plyr)
library(ggplot2)


# Globals -----------------------------------------------------------------


my.setups.brewer <- c("Standard"="#E41A1C",
                      "4 Input"="#377EB8", "6 Input"="#4DAF4A",
                      "10 Input"="#984EA3", "12 Input"="#FF7F00",
                      "2 Activated"="#FFFF33", "6 Activated"="#A65628",
                      "8 Activated"="#F781BF")
my.complexity.brewer <- c("2.2"="#E41A1C", "4.0"="#377EB8", "9.3"="#4DAF4A",
                          "12.9"="#984EA3")
my.setups.shapes <- c("Standard"=16, "4 Input"=4, "6 Input"=3, "10 Input"=11,
                      "12 Input"=8, "2 Activated"=15, "6 Activated"=17,
                      "8 Activated"=12)
my.complexity.shapes <- c("2.2"=16, "4.0"=15, "9.3"=17, "12.9"=18)
my.setups.lty <- c("Standard"=1, "4 Input"=2, "6 Input"=3, "10 Input"=4,
                      "12 Input"=5, "2 Activated"=2, "6 Activated"=3,
                      "8 Activated"=6)
my.complexity.lty <- c("2.2"=1, "4.0"=2, "9.3"=3, "12.9"=4)
my.milo.order <- sprintf("mtf_%d", 1:13)


# Utility Functions -------------------------------------------------------


probability_distributions <- function(my.df, my.sep, my.var, binw)
{
    return(ddply(my.df, my.sep, function(x) {
            x <- x[is.finite(x[[my.var]]),]
            tmp <- ggplot2:::bin(x[[my.var]], binwidth=binw, drop=TRUE,
                    range=c(min(x[[my.var]]) - binw,
                    max(x[[my.var]]) + binw))
            tmp[["y"]] <- tmp[["count"]] / sum(tmp[["count"]])
            for (my.name in my.sep) {
                label <- unique(x[[my.name]])
                tmp[[my.name]] <- rep(label, nrow(tmp))
            }
            tmp
        })
    )
}


# Journal Specific Output -------------------------------------------------


write_epjb_figure <- function(my.plot, my.path, my.title="", tall=1)
{
    # EPJB column width is 8.8 cm but in that size the standard font size is off
    my.w.cm <- 8.8 * 1.5
    # convert to standard inches
    my.w <- my.w.cm / 2.54
    my.h <- my.w * 3 / 4 * tall
    setEPS(reset=TRUE)
    postscript(file=paste(my.path, ".eps", sep=""), title=my.title,
               width=my.w, height=my.h)
    print(my.plot)
    dev.off()
    return()
}

write_normal <- function(my.plot, my.path, my.title="", tall=1)
{
    my.w <- 7
    my.h <- my.w * 3 / 4 * tall
    pdf(file=paste(my.path, ".pdf", sep=""), title=my.title,
        width=my.w, height=my.h)
    print(my.plot)
    dev.off()
    return()
}


# Plot Layouts ------------------------------------------------------------


layout_distribution <- function(my.plot,
                                my.xlab=expression(paste("X label")),
                                my.ylab=expression(paste("Y label")),
                                my.legend="Group",
                                my.a=1,
                                my.palette=my.setups.brewer,
                                my.lty=my.setups.lty,
                                my.l_sz=0.5)
{
    my.plot <- my.plot + geom_bar(stat="identity", position="identity",
                                  fill="transparent", alpha=my.a,
                                  show_guide=FALSE)
    my.plot <- my.plot + geom_line(alpha=my.a, size=my.l_sz)
    my.plot <- my.plot + scale_linetype_manual(name=my.legend, values=my.lty)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=my.palette)
    my.plot <- my.plot + scale_x_continuous(name=as.expression(bquote(.(my.xlab))))
    my.plot <- my.plot + scale_y_continuous(name=as.expression(bquote(.(my.ylab))))
    # prevent alpha values in the plot from reducing visibility of the legend
    my.plot <- my.plot + guides(colour=guide_legend(override.aes=list(alpha=1)))
    return(my.plot)
}

layout_scatter <- function(my.plot,
                           my.xlab=expression(paste("X label")),
                           my.ylab=expression(paste("Y label")),
                           my.legend="Group",
                           my.a=1,
                           my.palette=my.setups.brewer,
                           my.shapes=my.setups.shapes,
                           my.p_sz=1)
{
    my.plot <- my.plot + geom_point(alpha=my.a, size=my.p_sz)
    my.plot <- my.plot + scale_shape_manual(name=my.legend, values=my.shapes)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=my.palette)
    my.plot <- my.plot + scale_x_continuous(name=as.expression(bquote(.(my.xlab))))
    my.plot <- my.plot + scale_y_continuous(name=as.expression(bquote(.(my.ylab))))
    # prevent alpha values in the plot from reducing visibility of the legend
    my.plot <- my.plot + guides(colour=guide_legend(override.aes=list(alpha=1)))
    return(my.plot)
}

layout_tsp_with_error <- function(my.plot,
                                  my.xlab=expression(paste("Triad")),
                                  my.ylab=expression(paste("Z-Score")),
                                  my.legend="Group",
                                  my.a=1,
                                  my.palette=my.setups.brewer,
                                  my.shapes=my.setups.shapes,
                                  my.lty=my.setups.lty,
                                  my.p_sz=2,
                                  my.l_sz=0.5)
{
    my.plot <- my.plot + geom_point(alpha=my.a, size=my.p_sz)
    my.plot <- my.plot + geom_line(alpha=my.a, size=my.l_sz)
    my.plot <- my.plot + geom_errorbar(mapping=aes(linetype=NULL), alpha=my.a,
                                       width=0.5, size=my.l_sz, show_guide=FALSE)
    my.plot <- my.plot + scale_linetype_manual(name=my.legend, values=my.lty)
    my.plot <- my.plot + scale_shape_manual(name=my.legend, values=my.shapes)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=my.palette)
    my.plot <- my.plot + scale_x_discrete(limits=my.milo.order,
        name=as.expression(bquote(.(my.xlab))), labels=1:13)
    my.plot <- my.plot + scale_y_continuous(name=as.expression(bquote(.(my.ylab))))
    # prevent alpha values in the plot from reducing visibility of the legend
    my.plot <- my.plot + guides(colour=guide_legend(override.aes=list(alpha=1)))
    return(my.plot)
}


# Plotting ----------------------------------------------------------------


plot_all <- function(my.df, my.dest, my.title, my.legend, my.palette, my.shapes,
                     my.lty, my.sep=c("type", "setup"), my.facet="type ~ .",
                     my.write=write_normal, tall=2)
{
    use.facet <- TRUE
    t <- try(as.formula(my.facet), silent=TRUE)
    if (inherits(t, "try-error")) {
        use.facet <- FALSE
    }
    # distributions
    my.plot <- plot_flow_error(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "flow_error"),
             paste(my.title, "Flow Error"), tall=tall)

    my.plot <- plot_robustness(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "robustness"),
             paste(my.title, "Robustness"), tall=tall)

    my.plot <- plot_overlap(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "overlap"),
             paste(my.title, "Overlap"), tall=tall)

    my.plot <- plot_variance(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "pattern_variance"),
             paste(my.title, "Pattern Variance"), tall=tall)

    my.plot <- plot_scalar_complexity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "scalar_complexity"),
             paste(my.title, "Scalar Complexity"), tall=tall)

    my.plot <- plot_binary_complexity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "binary_complexity"),
             paste(my.title, "Binary Complexity"), tall=tall)
    
    my.plot <- plot_variance(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "pattern_variance"),
             paste(my.title, "Pattern Variance"), tall=tall)

    my.plot <- plot_binary_rank(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "binary_rank"),
             paste(my.title, "Binary Rank"), tall=tall)
    
    my.plot <- plot_pattern_rank(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "pattern_rank"),
             paste(my.title, "Pattern Rank"), tall=tall)
    
    my.plot <- plot_iteration(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "iteration"),
             paste(my.title, "Iteration"), tall=tall)

    my.plot <- plot_connectivity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "connectivity"),
             paste(my.title, "Connectivity"), tall=tall)

    my.plot <- plot_initial_connectivity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "initial_connectivity"),
             paste(my.title, "Initial Connectivity"), tall=tall)

    my.plot <- plot_spectral_modularity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "spectral_modularity"),
             paste(my.title, "Spectral Modularity"), tall=tall)

    my.plot <- plot_louvain_modularity(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "louvain_modularity"),
             paste(my.title, "Louvain Modularity"), tall=tall)

    my.plot <- plot_degree_correlation(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "degree_correlation"),
             paste(my.title, "Degree Correlation"), tall=tall)

    my.plot <- plot_average_shortest_path(my.df, my.legend, my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "average_shortest_path"),
             paste(my.title, "Average Shortest Path"), tall=tall)

    my.plot <- plot_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7",
                           my.palette, my.sep, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "zscore_triad_7"),
             paste(my.title, "Z-Score Triad 7"), tall=tall)
    # scatter plots
    my.plot <- plot_variance_vs_connectivity(my.df, my.legend, my.palette,
                                             my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "variance_vs_connectivity"),
             paste(my.title, "Pattern Variance vs Connectivity"), tall=tall)

    my.plot <- plot_binary_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                      my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "binary_vs_overlap"),
             paste(my.title, "Binary Complexity vs Overlap"), tall=tall)

    my.plot <- plot_binary_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                       my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "binary_vs_spectral"),
             paste(my.title, "Binary Complexity vs Modularity"), tall=tall)

    my.plot <- plot_scalar_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                      my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "scalar_vs_overlap"),
             paste(my.title, "Scalar Complexity vs Overlap"), tall=tall)

    my.plot <- plot_scalar_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                       my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "scalar_vs_spectral"),
             paste(my.title, "Scalar Complexity vs Modularity"), tall=tall)

    my.plot <- plot_scalar_vs_degree_corr(my.df, my.legend, my.palette, my.shapes,
                                          my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "scalar_vs_degree"),
             paste(my.title, "Scalar Complexity vs Degree Correlation"), tall=tall)

    my.plot <- plot_scalar_vs_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7",
                                     my.palette, my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "scalar_vs_zscore_triad_7"),
             paste(my.title, "Scalar Complexity vs Z-Score Triad 7"), tall=tall)
    
    my.plot <- plot_variance_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                         my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "variance_vs_spectral"),
             paste(my.title, "Pattern Variance vs Modularity"), tall=tall)

    my.plot <- plot_binary_rank_vs_spectral(my.df, my.legend, my.palette,
                                            my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "binary_rank_vs_spectral"),
             paste(my.title, "Binary Rank vs Modularity"), tall=tall)
    
    my.plot <- plot_pattern_rank_vs_spectral(my.df, my.legend, my.palette,
                                             my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "pattern_rank_vs_spectral"),
             paste(my.title, "Pattern Rank vs Modularity"), tall=tall)
    
    my.plot <- plot_spectral_vs_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7",
                                       my.palette, my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "spectral_vs_zscore_triad_7"),
             paste(my.title, "Modularity vs Z-Score Triad 7"), tall=tall)

    my.plot <- plot_spectral_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                        my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "spectral_vs_overlap"),
             paste(my.title, "Modularity vs Overlap"), tall=tall)

    my.plot <- plot_spectral_vs_degree_corr(my.df, my.legend, my.palette,
                                            my.shapes, my.a=1/5)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet)
    }
    my.write(my.plot, file.path(my.dest, "spectral_vs_degree"),
             paste(my.title, "Modularity vs Degree Correlation"), tall=tall)
#    # tsps
    tmp <- load_tsps(my.df, my.sep)
    my.plot <- plot_tsp(tmp, my.legend, my.palette, my.shapes, my.lty)
    if (use.facet) {
        my.plot <- my.plot + facet_grid(my.facet, scale="free_y")
    }
    my.write(my.plot, file.path(my.dest, "tsps"),
             paste(my.title, "Triad Significance Profiles"), tall=tall)
    return()
}

plot_publication <- function(my.dest, my.write=write_epjb_figure)
{
    my.legend <- "Setup"
    my.title <- "Parameter Setup, Phase I,"
    my.plot <- plot_scalar_complexity(my.setups.p1, my.legend, my.setups.brewer,
                                      my.sep=c("setup"), my.setups.lty, binw=0.2)
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "setups_phase_1_scalar_complexity"),
             paste(my.title, "Scalar Complexity"), tall=1)
    my.rob <- my.setups.p2[my.setups.p2$robustness > 0.6,]
    #     my.df <- subset(my.setups.p2, setup %in% c("6 Activated", "Standard",
    #             "2 Activated", "4 Input", "12 Input"))
    tmp <- subset(my.rob, setup %in% c("6 Activated", "Standard", "2 Activated"))
    #     tmp$setup <- factor(tmp$setup)
    tmp$setup <- factor(tmp$setup, levels=c("Standard", "2 Activated", "6 Activated"))
    my.title <- "Parameter Setup, Phase II,"
    # 
    my.plot <- plot_scalar_vs_spectral(tmp, my.legend, my.setups.brewer,
                                       my.setups.shapes, my.a=1, my.p_sz=1.5)
    my.plot <- my.plot + facet_grid("type ~ .")
#     my.plot <- my.plot + geom_smooth(mapping=aes(linetype=setup, group=setup),
#                                      method="lm", show_guide=FALSE)
#     my.plot <- my.plot + scale_linetype_manual(name=my.legend, values=my.setups.lty)
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "setups_phase_2_subset_scalar_vs_spectral"),
             paste(my.title, "Scalar Complexity vs Modularity"), tall=1.5)
    my.legend <- "Complexity"
    my.title <- "Artificial Patterns,"
    tmp <- subset(my.comp.p1, type == "Link Robust")
    my.total <- add_label(tmp, "phase", "Phase I")
    tmp <- subset(my.comp.p2, type == "Link Robust", robustness > 0.6)
    tmp <- add_label(tmp, "phase", "Phase II")
    my.total <- rbind(my.total, tmp)
    my.total$type <- factor(my.total$type)
    my.plot <- plot_spectral_modularity(my.total, my.legend, my.complexity.brewer,
                                        my.sep=c("setup", "phase"),
                                        my.complexity.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("phase ~ .")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_spectral_modularity"),
             paste(my.title, "Modularity"), tall=1)
    my.title <- "Artificial Patterns, Phase I,"
    tmp <- load_tsps(my.comp.p1, "setup")
    my.plot <- plot_tsp(tmp, my.legend, my.complexity.brewer,
                        my.complexity.shapes, my.complexity.lty)
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_phase_1_tsps"),
             paste(my.title, "Triad Significance Profiles"), tall=3/4)
    my.title <- "Artificial Patterns, Phase II,"
    tmp <- subset(my.comp.p2, setup %in% c("4.0", "9.3") & robustness > 0.6)
    tmp$setup <- factor(tmp$setup)
    tmp <- load_tsps(tmp, c("type", "setup"))
    my.plot <- plot_tsp(tmp, my.legend, my.complexity.brewer,
                        my.complexity.shapes, my.complexity.lty)
    my.plot <- my.plot + facet_grid("type ~ .", scale="free_y")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_phase_2_subset_tsps"),
             paste(my.title, "Triad Significance Profiles"), tall=1.5)
    # artificial patterns
    my.plot <- plot_patterns()
    my.write(my.plot, file.path(my.dest, "artificial_patterns"),
             paste(my.title, "Artificial Output Patterns"), tall=4/3)
}

#TODO
plot_appendix <- function(my.dest, my.write=write_epjb_figure)
{
    my.legend <- "Complexity"
    my.title <- "Artificial Patterns, Phase II"
    my.rob <- subset(my.comp.p2, robustness > 0.6)
    my.rob$setup <- factor(my.rob$setup)
    # modularity
    my.plot <- plot_spectral_modularity(my.rob, my.legend,
                                        my.complexity.brewer,
                                        my.sep=c("type", "setup"),
                                        my.complexity.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("type ~ .")
#     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_phase_2_spectral_modularity"),
             paste(my.title, "Modularity"), tall=1.5)
    # density
    my.plot <- plot_connectivity(my.rob, my.legend,
                                        my.complexity.brewer,
                                        my.sep=c("type", "setup"),
                                        my.complexity.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("type ~ .")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_phase_2_connectivity"),
             paste(my.title, "Connectivity"), tall=1.5)
    # tsps
    tmp <- load_tsps(my.rob, c("type", "setup"))
    my.plot <- plot_tsp(tmp, my.legend, my.complexity.brewer,
                        my.complexity.shapes, my.complexity.lty)
    my.plot <- my.plot + facet_grid("type ~ .", scale="free_y")
#     print(my.plot)
    my.write(my.plot, file.path(my.dest, "artificial_phase_2_tsps"),
             paste(my.title, "Triad Significance Profiles"), tall=1.5)
    my.legend <- "Setup"
    my.title <- "Parameter Setups, Phase II"
    my.rob <- subset(my.setups.p2, robustness > 0.6)
    my.rob$setup <- factor(my.rob$setup)
    # binary rank
    my.plot <- plot_binary_rank(my.rob, my.legend,
                                 my.setups.brewer,
                                 my.sep=c("type", "setup"),
                                 my.setups.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("type ~ .")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "setups_phase_2_binary_rank"),
             paste(my.title, "Binary Rank"), tall=1.5)
    # modularity
    my.plot <- plot_spectral_modularity(my.rob, my.legend,
                            my.setups.brewer,
                            my.sep=c("type", "setup"),
                            my.setups.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("type ~ .")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "setups_phase_2_spectral_modularity"),
             paste(my.title, "Spectral Modularity"), tall=1.5)
    # mean overlap
    my.plot <- plot_overlap(my.rob, my.legend,
                                 my.setups.brewer,
                                 my.sep=c("type", "setup"),
                                 my.setups.lty, binw=0.01)
    my.plot <- my.plot + facet_grid("type ~ .")
    #     print(my.plot)
    my.write(my.plot, file.path(my.dest, "setups_phase_2_mean_overlap"),
             paste(my.title, "Mean Overlap"), tall=1.5)
}


# Distributions of Single Attributes --------------------------------------


plot_flow_error <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                            binw=10^-4, my.a=1, lim=0.0075, thresh=0.007)
{
    cat("flow error\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    cat(paste("\tright limit =", lim, "\n"))
    cat(paste("\tthreshold =", thresh, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "flow_error", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Flow Error ", epsilon)),
                                   expression(paste(PMF(epsilon))),
                                   my.legend, my.a, my.palette, my.lty)
    my.plot <- my.plot + coord_cartesian(xlim=c(0, lim))
    # mark threshold in plot
    my.plot <- my.plot + geom_vline(xintercept=thresh, linetype=1, size=0.3)
    # make dataframe for label
    if (length(my.sep) > 1) {
        my.text <- ddply(tmp, "type", colwise(max, "y"))
        my.text$x <- rep(thresh, length(levels(tmp$type)))
    }
    else {
        my.text <- data.frame(x=thresh, y=max(tmp$y))
    }
    my.text$y <- my.text$y + my.text$y / 20
    my.plot <- my.plot + geom_text(data=my.text, mapping=aes(colour=NULL,
                                   linetype=NULL), colour="black", hjust=1.1,
            label="paste('Threshold ', italic(h))", size=4, parse=TRUE)
    return(my.plot)
}

plot_robustness <- function(my.df, my.legend, my.palette, my.sep, my.lty, binw=0.01, my.a=1)
{
    cat("robustness\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "robustness", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Robustness ", rho)),
                                   expression(paste(PMF(rho))), my.legend, my.a,
                                   my.palette, my.lty)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    return(my.plot)
}

plot_overlap <- function(my.df, my.legend, my.palette, my.sep, my.lty, binw=0.01, my.a=1)
{
    cat("overlap\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "mean_overlap", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Overlap ", O)),
                                   expression(paste(PMF(O))), my.legend, my.a,
                                   my.palette, my.lty)
#     my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    return(my.plot)
}

plot_variance <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                          binw=0.01, my.a=1)
{
    cat("variance\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "pattern_variance", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Pattern Variance ", sigma^2[Q])),
                                 expression(paste(PMF(sigma^2[Q]))), my.legend,
                                   my.a, my.palette, my.lty)
    return(my.plot)
}

plot_binary_rank <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                             binw=0.1, my.a=1)
{
    cat("binary rank\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "binary_rank", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Binary Rank ", R[b])),
                                   expression(paste(PMF(R[b]))), my.legend,
                                   my.a, my.palette, my.lty)
    return(my.plot)
}

plot_pattern_rank <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                              binw=0.1, my.a=1)
{
    cat("pattern rank\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "pattern_rank", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Pattern Rank ", R[Q])),
                                   expression(paste(PMF(R[Q]))), my.legend,
                                   my.a, my.palette, my.lty)
    return(my.plot)
}

plot_scalar_complexity <- function(my.df, my.legend, my.palette, my.sep,
                                   my.lty, binw=0.01, my.a=1)
{
    cat("scalar complexity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "scalar_complexity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Scalar Complexity ", italic(C))),
                                 expression(paste(PMF(italic(C)))), my.legend,
                                   my.a, my.palette, my.lty)
    return(my.plot)
}

plot_binary_complexity <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                                   binw=0.05, my.a=1)
{
    cat("binary complexity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "binary_complexity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Binary Complexity ", italic(C[b]))),
                                 expression(paste(PMF(italic(C[b])))),
                                   my.legend, my.a, my.palette, my.lty)
    return(my.plot)
}

plot_iteration <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                           binw=10, my.a=1, lim=10^4)
{
    cat("iteration\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    cat(paste("\tright limit =", lim, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "iteration", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Iteration ", italic(i))),
                                 expression(paste(PMF(italic(i)))),
                                   my.legend, my.a, my.palette, my.lty)
    my.plot <- my.plot + coord_cartesian(xlim=c(0, lim))
    return(my.plot)
}

plot_connectivity <- function(my.df, my.legend, my.palette, my.sep, my.lty,
                              binw=0.001, my.a=1)
{
    cat("connectivity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "density", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Connectivity ", gamma)),
                                 expression(paste(PMF(gamma))), my.legend,
                                   my.a, my.palette, my.lty)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    return(my.plot)
}

plot_initial_connectivity <- function(my.df, my.legend, my.palette, my.sep,
                                      my.lty, binw=0.01, my.a=1)
{
    cat("initial connectivity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep,
            "initial_connectivity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Initial Connectivity ", gamma[init])),
                                   expression(paste(PMF(gamma[init]))),
                                   my.legend, my.a, my.palette, my.lty)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    return(my.plot)
}

plot_spectral_modularity <- function(my.df, my.legend, my.palette, my.sep,
                                     my.lty, binw=0.01, my.a=1)
{
    cat("spectral modularity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep,
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Modularity ", italic(Q))),
                                   expression(paste(PMF(italic(Q)))),
                                   my.legend, my.a, my.palette, my.lty)
    return(my.plot)
}

plot_louvain_modularity <- function(my.df, my.legend, my.palette, my.sep,
                                    my.lty, binw=0.01, my.a=1)
{
    cat("louvain modularity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep,
            "louvain_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Modularity ", italic(Q))),
                                   expression(paste(PMF(italic(Q)))),
                                   my.legend, my.a, my.palette, my.lty)
    return(my.plot)
}

plot_degree_correlation <- function(my.df, my.legend, my.palette, my.sep,
                                    my.lty, binw=0.01, my.a=1)
{
    cat("degree correlation\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "degree_correlation", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Degree Correlation ", italic(r))),
                                   expression(paste(PMF(italic(r)))),
                                   my.legend, my.a, my.palette, my.lty)
#    my.plot <- my.plot + coord_cartesian(xlim=c(-1, 1))
    return(my.plot)
}

plot_average_shortest_path <- function(my.df, my.legend, my.palette, my.sep,
                                       my.lty, binw=0.01, my.a=1)
{
    cat("shortest path\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, "shortest_paths", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Average Shortest Path ", italic(p))),
                                   expression(paste(PMF(italic(p)))),
                                   my.legend, my.a, my.palette, my.lty)
    return(my.plot)
}

plot_zscore <- function(my.df, my.legend, mtf, label, my.palette, my.sep,
                        my.lty, binw=0.2, my.a=1)
{
    cat(paste(bquote(.(mtf)), "\n"))
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, my.sep, mtf, binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, label,
                                   as.expression(bquote(PMF(.(label)))),
                                   my.legend, my.a, my.palette, my.lty)
    return(my.plot)
}


# Scatter Plots -----------------------------------------------------------


plot_comparison <- function(my.df, my.attr, my.legend, my.palette, my.shapes, my.a=1, axis.sz=6)
{
    cat("plotting comparison\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=phase_1, y=phase_2, colour=setup))
    my.plot <- layout_scatter(my.plot, paste(my.attr, "Phase I"),
                              paste(my.attr, "Phase II"),
                              my.legend, my.a, my.palette, my.shapes)
    my.plot <- my.plot + geom_abline(intercept=0, slope=1, colour="red",
                                     show_guide=FALSE)
    return(my.plot)
}

plot_variance_vs_connectivity <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("variance vs density\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=pattern_variance, y=density,
                                 colour=setup, linetype=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Pattern Variance ", sigma^2)),
                              expression(paste("Connectivity ", gamma)),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_binary_vs_overlap <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("binary vs overlap\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=binary_complexity, y=mean_overlap, colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Binary Complexity ", italic(C[b]))),
                              expression(paste("Overlap ", italic(O))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_binary_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("binary vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=binary_complexity, y=spectral_modularity,
                                 colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Binary Complexity ", italic(C[b]))),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_scalar_vs_overlap <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("scalar vs overlap\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=scalar_complexity, y=mean_overlap,
                                 colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Scalar Complexity ", italic(C))),
                              expression(paste("Overlap ", italic(O))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_scalar_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes, my.a=1, my.p_sz=1.5)
{
    cat("scalar vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=scalar_complexity, y=spectral_modularity,
                                     colour=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Scalar Complexity ", italic(C))),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}

plot_scalar_vs_degree_corr <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("scalar vs degree correlation\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=scalar_complexity, y=degree_correlation,
                                 colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Scalar Complexity ", italic(C))),
                              expression(paste("Degree Correlation ", italic(r))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_scalar_vs_zscore <- function(my.df, my.legend, mtf, label, my.palette, my.shapes, my.a=1)
{
    cat(paste("scalar vs", mtf, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes_string(x="scalar_complexity", y=mtf,
                                        colour="setup"))
    my.plot <- layout_scatter(my.plot,
                              label,
                              expression(paste("Scalar Complexity ", italic(C))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_variance_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes, my.a=1, my.p_sz=1.5)
{
    cat("variance vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=pattern_variance, y=spectral_modularity,
                                 colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Pattern Variance ", sigma^2)),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}

plot_variance_vs_overlap <- function(my.df, my.legend, my.palette, my.shapes, my.a=1, my.p_sz=1.5)
{
    cat("variance vs overlap\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=pattern_variance, y=mean_overlap,
                                 colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Pattern Variance ", sigma^2)),
                              expression(paste("Overlap ", italic(O))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}


plot_binary_rank_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes,
                                         my.a=1, my.p_sz=1.5)
{
    cat("binary rank vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=binary_rank, y=spectral_modularity,
                                 colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Binary Rank ", R[b])),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}

plot_binary_rank_vs_overlap <- function(my.df, my.legend, my.palette, my.shapes,
                                         my.a=1, my.p_sz=1.5)
{
    cat("binary rank vs overlap\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=binary_rank, y=mean_overlap,
                                 colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Binary Rank ", R[b])),
                              expression(paste("Overlap ", italic(O))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}

plot_pattern_rank_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes,
                                  my.a=1, my.p_sz=1.5)
{
    cat("pattern rank vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=pattern_rank, y=spectral_modularity,
                                 colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Pattern Rank ", R[Q])),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    return(my.plot)
}

plot_spectral_vs_zscore <- function(my.df, my.legend, mtf, label, my.palette, my.shapes, my.a=1)
{
    cat(paste("modularity vs", mtf, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes_string(x="spectral_modularity", y=mtf,
                                 colour="setup"))
    my.plot <- layout_scatter(my.plot,
                              label,
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_spectral_vs_overlap <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("modularity vs overlap\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=spectral_modularity, y=mean_overlap,
                                 colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Modularity ", italic(Q))),
                              expression(paste("Overlap ", italic(O))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_spectral_vs_degree_corr <- function(my.df, my.legend, my.palette, my.shapes, my.a=1)
{
    cat("modularity vs degree correlation\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=spectral_modularity, y=degree_correlation,
                                 colour=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Modularity ", italic(Q))),
                              expression(paste("Degree Correlation ", italic(r))),
                              my.legend, my.a, my.palette, my.shapes)
    return(my.plot)
}

plot_overlap_modularity_boundaries <- function(my.legend="Setup", my.a=1)
{
    final <- subset(my.setups.final, robustness > 0.6)
    # overall data frame
    tmp <- subset(final, type=="Link Robust")
    res <- subset(tmp, setup %in% c("6 Activated", "Standard", "2 Activated"))
    tmp <- subset(final, type=="Node Robust")
    tmp <- subset(tmp, setup %in% c("6 Activated", "2 Activated", "4 Input"))
    res <- rbind(res, tmp)
    tmp <- subset(final, type=="Noise Robust")
    tmp <- subset(tmp, setup %in% c("6 Activated", "Standard", "4 Input"))
    res <- rbind(res, tmp)
    res$setup <- factor(res$setup)
    # data frame for linear fit
    tmp <- subset(final, type=="Link Robust")
    fit <- subset(tmp, setup == "Standard")
    tmp <- subset(final, type == "Node Robust")
    tmp <- subset(tmp, setup == "6 Activated")
    fit <- rbind(fit, tmp)
    tmp <- subset(final, type=="Noise Robust")
    tmp <- subset(tmp, setup == "6 Activated")
    fit <- rbind(fit, tmp)
    fit$setup <- factor(fit$setup)
    # consistent colours
    colours <- RColorBrewer::brewer.pal(8, my.brewer)
    my.plot <- ggplot(res)
    my.plot <- my.plot + geom_point(mapping=aes(x=spectral_modularity,
                                    y=mean_middle_overlap, colour=setup, shape=setup),
                                    size=1.5, alpha=my.a)
    my.plot <- my.plot + scale_shape(name=my.legend)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=c(colours[1],
                                             colours[2], colours[4], colours[5]))
    my.plot <- my.plot + facet_grid("type ~ .")
    my.plot <- my.plot + scale_y_continuous(
        name=expression(paste("Overlap ", italic(O))))
    my.plot <- my.plot + scale_x_continuous(
        name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + geom_smooth(data=fit, mapping=aes(x=spectral_modularity,
            y=mean_middle_overlap, group=setup), colour="black", method="lm")
#     my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot
}


# TSPs --------------------------------------------------------------------


plot_tsp <- function(my.df, my.legend, my.palette, my.shapes, my.lty, my.a=1)
{
    cat("tsps\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=triad, y=zscore, ymax=zscore + error,
                               ymin=zscore - error, colour=setup, group=setup,
                               linetype=setup, shape=setup))
    my.plot <- layout_tsp_with_error(my.plot,
                                     my.xlab=expression(paste("Triad")),
                                     my.ylab=expression(paste("Z-Score")),
                                     my.legend,
                                     my.a,
                                     my.palette,
                                     my.shapes,
                                     my.lty)
    my.plot
}


# Output Patterns ---------------------------------------------------------


plot_patterns <- function()
{
    my.df <- read.csv("complexity/output_patterns.csv")
    my.plot <- ggplot(my.df, aes(x=j, y=i, fill=value))
    my.plot <- my.plot + theme_bw()
    my.plot <- my.plot + geom_tile(colour="black")
    my.plot <- my.plot + facet_wrap(~ setup, ncol=2)
    my.plot <- my.plot + scale_fill_gradient(name="", low="white",
            high="black", limits=c(0, 1))
    my.plot <- my.plot + theme(axis.ticks=element_blank(),
            panel.background=element_blank(), plot.background=element_blank(),
            panel.grid.minor=element_blank(), panel.grid.major=element_blank()
            )
    my.plot <- my.plot + scale_y_reverse(name=
            expression(paste("Output Node ", italic(i))), breaks=1:max(my.df$i))
    my.plot <- my.plot + scale_x_continuous(name=
            expression(paste("Input Node ", italic(j))), breaks=1:max(my.df$j))
    return(my.plot)
}
