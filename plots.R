
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


# Utility Functions -------------------------------------------------------


probability_distributions <- function(my.df, my.sep, my.var, binw)
{
    return(ddply(my.df, my.sep, function(x) {
            tmp <- ggplot2:::bin(x[[my.var]], binwidth=binw, drop=TRUE,
                    range=c(min(x[[my.var]], na.rm=TRUE) - binw,
                    max(x[[my.var]], na.rm=TRUE) + binw))
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


write_epjb_figure <- function(my.plot, my.path, my.title="", tall=FALSE)
{
    # EPJB column width is 8.8 cm but in that size the standard font size is off
    my.w.cm <- 8.8 * 1.5
    my.w <- my.w.cm / 2.54
    if (tall) {
        my.h <- my.w * 1.5
    }
    else {
        my.h <- my.w
    }
    setEPS(reset=TRUE)
    postscript(file=paste(my.path, ".eps", sep=""), title=my.title, width=my.w, height=my.h)
    print(my.plot)
    dev.off()
    return()
}

write_normal <- function(my.plot, my.path, my.title="", tall=FALSE)
{
    my.w <- 7
    if (tall) {
        my.h <- my.w * 1.5
    }
    else {
        my.h <- my.w
    }
    pdf(file=paste(my.path, ".pdf", sep=""), title=my.title, width=my.w, height=my.h)
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
                                my.l_sz=0.5)
{
    my.plot <- my.plot + geom_bar(stat="identity", position="identity",
                                  fill="transparent", alpha=my.a,
                                  show_guide=FALSE)
    my.plot <- my.plot + geom_line(alpha=my.a, size=my.l_sz)
    my.plot <- my.plot + scale_linetype(name=my.legend)
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
                                  my.p_sz=2,
                                  my.l_sz=0.5)
{
    my.plot <- my.plot + geom_point(alpha=my.a, size=my.p_sz)
    my.plot <- my.plot + geom_line(alpha=my.a, size=my.l_sz)
    my.plot <- my.plot + geom_errorbar(alpha=my.a, size=my.l_sz, width=0.4,
                                       show_guide=FALSE)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_shape_manual(name=my.legend, values=my.shapes)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=my.palette)
    my.plot <- my.plot + scale_x_continuous(
        name=as.expression(bquote(.(my.xlab))), breaks=1:13)
    # prevent alpha values in the plot from reducing visibility of the legend
    my.plot <- my.plot + guides(colour=guide_legend(override.aes=list(alpha=1)))
    return(my.plot)
}


# Plotting ----------------------------------------------------------------


plot_all <- function(my.df, my.dest, my.title, my.legend, my.palette, my.shapes,
                     my.write=write_normal)
{
    # distributions
    my.plot <- plot_flow_error(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "flow_error"),
             paste(my.title, "Flow Error"), tall=TRUE)

    my.plot <- plot_robustness(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "robustness"),
             paste(my.title, "Robustness"), tall=TRUE)

    my.plot <- plot_overlap(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "overlap"),
             paste(my.title, "Overlap"), tall=TRUE)

    my.plot <- plot_variance(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "pattern_variance"),
             paste(my.title, "Pattern Variance"), tall=TRUE)

    my.plot <- plot_scalar_complexity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "scalar_complexity"),
             paste(my.title, "Scalar Complexity"), tall=TRUE)

    my.plot <- plot_binary_complexity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "binary_complexity"),
             paste(my.title, "Binary Complexity"), tall=TRUE)
    
    my.plot <- plot_variance(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "pattern_variance"),
             paste(my.title, "Pattern Variance"), tall=TRUE)

    my.plot <- plot_binary_rank(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "binary_rank"),
             paste(my.title, "Binary Rank"), tall=TRUE)
    
    my.plot <- plot_pattern_rank(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "pattern_rank"),
             paste(my.title, "Pattern Rank"), tall=TRUE)
    
    my.plot <- plot_iteration(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "iteration"),
             paste(my.title, "Iteration"), tall=TRUE)

    my.plot <- plot_connectivity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "connectivity"),
             paste(my.title, "Connectivity"), tall=TRUE)

    my.plot <- plot_initial_connectivity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "initial_connectivity"),
             paste(my.title, "Initial Connectivity"), tall=TRUE)

    my.plot <- plot_spectral_modularity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "spectral_modularity"),
             paste(my.title, "Spectral Modularity"), tall=TRUE)

    my.plot <- plot_louvain_modularity(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "louvain_modularity"),
             paste(my.title, "Louvain Modularity"), tall=TRUE)

    my.plot <- plot_degree_correlation(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "degree_correlation"),
             paste(my.title, "Degree Correlation"), tall=TRUE)

    my.plot <- plot_average_shortest_path(my.df, my.legend, my.palette)
    my.write(my.plot, file.path(my.dest, "average_shortest_path"),
             paste(my.title, "Average Shortest Path"), tall=TRUE)

    my.plot <- plot_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7", my.palette)
    my.write(my.plot, file.path(my.dest, "zscore_triad_7"),
             paste(my.title, "Z-Score Triad 7"), tall=TRUE)
    # scatter plots
    my.plot <- plot_variance_vs_connectivity(my.df, my.legend, my.palette,
                                             my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "variance_vs_connectivity"),
             paste(my.title, "Pattern Variance vs Connectivity"), tall=TRUE)

    my.plot <- plot_binary_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                      my.a=1/5)
    my.write(my.plot, file.path(my.dest, "binary_vs_overlap"),
             paste(my.title, "Binary Complexity vs Overlap"), tall=TRUE)

    my.plot <- plot_binary_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                       my.a=1/5)
    my.write(my.plot, file.path(my.dest, "binary_vs_spectral"),
             paste(my.title, "Binary Complexity vs Modularity"), tall=TRUE)

    my.plot <- plot_scalar_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                      my.a=1/5)
    my.write(my.plot, file.path(my.dest, "scalar_vs_overlap"),
             paste(my.title, "Scalar Complexity vs Overlap"), tall=TRUE)

    my.plot <- plot_scalar_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                       my.a=1/5)
    my.write(my.plot, file.path(my.dest, "scalar_vs_spectral"),
             paste(my.title, "Scalar Complexity vs Modularity"), tall=TRUE)

    my.plot <- plot_scalar_vs_degree_corr(my.df, my.legend, my.palette, my.shapes,
                                          my.a=1/5)
    my.write(my.plot, file.path(my.dest, "scalar_vs_degree"),
             paste(my.title, "Scalar Complexity vs Degree Correlation"), tall=TRUE)

    my.plot <- plot_scalar_vs_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7",
                                     my.palette, my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "scalar_vs_zscore_triad_7"),
             paste(my.title, "Scalar Complexity vs Z-Score Triad 7"), tall=TRUE)
    
    my.plot <- plot_variance_vs_spectral(my.df, my.legend, my.palette, my.shapes,
                                         my.a=1/5)
    my.write(my.plot, file.path(my.dest, "variance_vs_spectral"),
             paste(my.title, "Pattern Variance vs Modularity"), tall=TRUE)

    my.plot <- plot_binary_rank_vs_spectral(my.df, my.legend, my.palette,
                                            my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "binary_rank_vs_spectral"),
             paste(my.title, "Binary Rank vs Modularity"), tall=TRUE)
    
    my.plot <- plot_pattern_rank_vs_spectral(my.df, my.legend, my.palette,
                                             my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "pattern_rank_vs_spectral"),
             paste(my.title, "Pattern Rank vs Modularity"), tall=TRUE)
    
    my.plot <- plot_spectral_vs_zscore(my.df, my.legend, "mtf_7", "Z-Score Triad 7",
                                       my.palette, my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "spectral_vs_zscore_triad_7"),
             paste(my.title, "Modularity vs Z-Score Triad 7"), tall=TRUE)

    my.plot <- plot_spectral_vs_overlap(my.df, my.legend, my.palette, my.shapes,
                                        my.a=1/5)
    my.write(my.plot, file.path(my.dest, "spectral_vs_overlap"),
             paste(my.title, "Modularity vs Overlap"), tall=TRUE)

    my.plot <- plot_spectral_vs_degree_corr(my.df, my.legend, my.palette,
                                            my.shapes, my.a=1/5)
    my.write(my.plot, file.path(my.dest, "spectral_vs_degree"),
             paste(my.title, "Modularity vs Degree Correlation"), tall=TRUE)

#    # tsps
#    pdf(file.path(my.dest, "tsps.pdf"),
#        title=paste(my.title, "TSPs"))
#    print(plot_tsps(my.setups.flow.tsp, my.legend))
#    dev.off()

    return()
}

plot_tsps_subset <- function(my.df, my.legend, a=1)
{
    print("tsps")
    print(paste("alpha =", a))
    # want to keep same colours as in all other plots
    colours <- RColorBrewer::brewer.pal(8, my.brewer)
    my.df <- subset(my.df, setup %in% c("4.0", "9.3"))
    my.df$setup <- factor(my.df$setup)
    tmp <- ddply(my.df, c("setup", "type", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE), uncrt=sd(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
                                   alpha=a)
    my.plot <- my.plot + geom_errorbar(mapping=aes(ymin=zscore - uncrt,
                                                   ymax=zscore + uncrt, colour=setup), alpha=a, legend=F)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_manual(name=my.legend, values=colours[2:3])
    my.plot <- my.plot + facet_grid("type ~ .", scales="free_y")
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")
    my.plot
}

plot_tsps_no_error <- function(my.df, my.legend, a=1)
{
    print("tsps")
    print(paste("alpha =", a))
    tmp <- ddply(my.df, c("setup", "type", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
            alpha=a)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_brewer(name=my.legend, palette=my.brewer)
    my.plot <- my.plot + facet_grid("type ~ .", scales="free_y")
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")
    my.plot
}


# Distributions of Single Attributes --------------------------------------


plot_flow_error <- function(my.df, my.legend, my.palette, binw=10^-4, my.a=1, lim=0.0075, thresh=0.007)
{
    cat("flow error\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    cat(paste("\tright limit =", lim, "\n"))
    cat(paste("\tthreshold =", thresh, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "flow_error", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Flow Error ", epsilon)),
                                   expression(paste(PMF(epsilon))),
                                   my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    my.plot <- my.plot + coord_cartesian(xlim=c(0, lim))
    # mark threshold in plot
    my.plot <- my.plot + geom_vline(xintercept=thresh, linetype=1, size=0.3)
    # make dataframe for label
    my.text <- ddply(tmp, "type", colwise(max, "y"))
    my.text$x <- rep(thresh, length(levels(tmp$type)))
    my.plot <- my.plot + geom_text(data=my.text, mapping=aes(colour=NULL,
                                   linetype=NULL), colour="black", hjust=1.1,
            label="paste('Threshold ', italic(h))", size=4, parse=TRUE)
    return(my.plot)
}

plot_robustness <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("robustness\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "robustness", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Robustness ", rho)),
                                   expression(paste(PMF(rho))), my.legend, my.a,
                                   my.palette)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_overlap <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("overlap\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "mean_overlap", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Overlap ", O)),
                                   expression(paste(PMF(O))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
#     my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    return(my.plot)
}

plot_variance <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("variance\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "pattern_variance", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Pattern Variance ", sigma^2[Q])),
                                 expression(paste(PMF(sigma^2[Q]))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_binary_rank <- function(my.df, my.legend, my.palette, binw=0.1, my.a=1)
{
    cat("binary rank\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "binary_rank", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Binary Rank ", R[b])),
                                   expression(paste(PMF(R[b]))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_pattern_rank <- function(my.df, my.legend, my.palette, binw=0.1, my.a=1)
{
    cat("binary rank\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "pattern_rank", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Pattern Rank ", R[Q])),
                                   expression(paste(PMF(R[Q]))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_scalar_complexity <- function(my.df, my.legend, my.palette, binw=0.005, my.a=1)
{
    cat("scalar complexity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "scalar_complexity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Scalar Complexity ", italic(C))),
                                 expression(paste(PMF(italic(C)))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_binary_complexity <- function(my.df, my.legend, my.palette, binw=0.05, my.a=1)
{
    cat("binary complexity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "binary_complexity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Binary Complexity ", italic(C[b]))),
                                 expression(paste(PMF(italic(C[b])))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_iteration <- function(my.df, my.legend, my.palette, binw=10, my.a=1, lim=10^4)
{
    cat("iteration\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    cat(paste("\tright limit =", lim, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "iteration", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Iteration ", italic(i))),
                                 expression(paste(PMF(italic(i)))), my.legend, my.a, my.palette)
    my.plot <- my.plot + coord_cartesian(xlim=c(0, lim))
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_connectivity <- function(my.df, my.legend, my.palette, binw=0.001, my.a=1)
{
    cat("connectivity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "density", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, expression(paste("Connectivity ", gamma)),
                                 expression(paste(PMF(gamma))), my.legend, my.a, my.palette)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_initial_connectivity <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("initial connectivity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"),
            "initial_connectivity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Initial Connectivity ", gamma[init])),
                                   expression(paste(PMF(gamma[init]))), my.legend, my.a, my.palette)
#    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1))
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_spectral_modularity <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("spectral modularity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Modularity ", italic(Q))),
                                   expression(paste(PMF(italic(Q)))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_louvain_modularity <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("louvain modularity\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"),
            "louvain_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Modularity ", italic(Q))),
                                   expression(paste(PMF(italic(Q)))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_degree_correlation <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("degree correlation\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "degree_correlation", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Degree Correlation ", italic(r))),
                                   expression(paste(PMF(italic(r)))), my.legend, my.a, my.palette)
#    my.plot <- my.plot + coord_cartesian(xlim=c(-1, 1))
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_average_shortest_path <- function(my.df, my.legend, my.palette, binw=0.01, my.a=1)
{
    cat("shortest path\n")
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), "shortest_paths", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot,
                                   expression(paste("Average Shortest Path ", italic(p))),
                                   expression(paste(PMF(italic(p)))), my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_zscore <- function(my.df, my.legend, mtf, label, my.palette, binw=0.2, my.a=1)
{
    cat(paste(bquote(.(mtf)), "\n"))
    cat(paste("\tbinwidth =", binw, "\n"))
    cat(paste("\talpha =", my.a, "\n"))
    tmp <- probability_distributions(my.df, c("setup", "type"), mtf, binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y, colour=setup, linetype=setup))
    my.plot <- layout_distribution(my.plot, label,
                                   as.expression(bquote(PMF(.(label)))),
                                   my.legend, my.a, my.palette)
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
    return(my.plot)
}

plot_scalar_vs_spectral <- function(my.df, my.legend, my.palette, my.shapes, my.a=1, my.p_sz=1.5)
{
    cat("scalar vs modularity\n")
    cat(paste("\talpha =", my.a, "\n"))
    my.plot <- ggplot(my.df, aes(x=scalar_complexity, y=spectral_modularity,
                                     colour=setup, linetype=setup, shape=setup))
    my.plot <- layout_scatter(my.plot,
                              expression(paste("Scalar Complexity ", italic(C))),
                              expression(paste("Modularity ", italic(Q))),
                              my.legend, my.a, my.palette, my.shapes, my.p_sz=my.p_sz)
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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
    my.plot <- my.plot + facet_grid("type ~ .")
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


# Output Patterns ---------------------------------------------------------


plot_patterns <- function()
{
    my.df <- read.csv("complexity/output_patterns.csv")
    my.plot <- ggplot(my.df, aes(x=j, y=i, fill=value))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_tile(colour="black")
    my.plot <- my.plot + facet_wrap(~ setup, ncol=2)
    my.plot <- my.plot + scale_fill_gradient(name="", low="white",
            high="black", limits=c(0, 1))
    my.plot <- my.plot + opts(axis.ticks=theme_blank(),
            panel.background=theme_blank(), plot.background=theme_blank(),
            panel.grid.minor=theme_blank(), panel.grid.major=theme_blank()
            )
    my.plot <- my.plot + scale_y_reverse(name=
            expression(paste("Output Node ", italic(i))), breaks=1:8)
    my.plot <- my.plot + scale_x_continuous(name=
            expression(paste("Input Node ", italic(j))), breaks=1:8)
    pdf(file.path("..", "Figures", "complexity", "output_patterns.pdf"),
            title="Artificial Output Patterns")
    cat(my.plot)
    dev.off()
}


# Statistics --------------------------------------------------------------


shifts <- function(my.df, my.var)
{
    my.names <- levels(my.df$setup)
    for (i in 1:(length(my.names) - 1)) {
        j <- i + 1
        cat(paste(my.names[i], my.names[j]))
        cat(t.test(subset(my.df, setup==my.names[i])[[my.var]],
                subset(my.df, setup==my.names[j])[[my.var]]))
    }
}


# Publication Plots ------------------------------------------------------------


plot_publication <- function(my.path, my.write=write_epjb_figure)
{
    my.df <- my.setups.p2[my.setups.p2$robustness > 0.6,]
#     my.df <- subset(my.setups.p2, setup %in% c("6 Activated", "Standard",
#             "2 Activated", "4 Input", "12 Input"))
    tmp <- subset(my.df, setup %in% c("6 Activated", "Standard", "2 Activated"))
    tmp$setup <- factor(tmp$setup)
    tmp$setup <- factor(tmp$setup, levels=c("Standard", "2 Activated", "6 Activated"))
    my.palette <- my.setups.brewer
    my.shapes <- my.setups.shapes
    my.legend <- "Parameter Setup"
    my.title <- "Parameter Setup, Phase II,"
    # 
    my.plot <- plot_scalar_vs_spectral(tmp, my.legend, my.palette, my.shapes, my.p_sz=3)
    my.plot <- my.plot + geom_smooth(method="lm") + scale_linetype(name=my.legend)
    print(my.plot)
#     my.write(my.plot, file.path(my.dest, "scalar_vs_spectral"),
#              paste(my.title, "Scalar Complexity vs Modularity"), tall=TRUE)
}

plot_figure_overlap <- function(my.legend="Setup", my.a=1)
{
    my.setups.final <- subset(my.setups.final, robustness > 0.6)
    tmp <- subset(my.setups.final, type=="Link Robust")
    res <- subset(tmp, setup == "Standard")
    tmp <- subset(my.setups.final, type == "Node Robust")
    tmp <- subset(tmp, setup == "6 Activated")
    res <- rbind(res, tmp)
    tmp <- subset(my.setups.final, type=="Noise Robust")
    tmp <- subset(tmp, setup == "6 Activated")
    res <- rbind(res, tmp)
    res$setup <- factor(res$setup)
    my.plot <- ggplot(res, aes(x=spectral_modularity, y=mean_middle_overlap, colour=type, shape=setup))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_point(size=2, alpha=my.a)
    my.plot <- my.plot + scale_colour_brewer(name="Type", palette=my.brewer)
    my.plot <- my.plot + scale_shape(name=my.legend)
#     my.plot <- my.plot + facet_grid("type ~ .")
    my.plot <- my.plot + scale_y_continuous(
        name=expression(paste("Overlap ", italic(o))))
    my.plot <- my.plot + scale_x_continuous(
        name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + geom_smooth(mapping=aes(colour=type, group=type,
                                                 linetype=type), method="lm")
    my.plot <- my.plot + scale_linetype(name="Type")
    my.plot
}
    
plot_additional <- function()
{
    my.dest <- file.path("..", "Figures")
    # plot everything for complexity after phase 1
#    my.sub <- file.path(my.dest, "complexity", "at_threshold")
#    my.title <- "Artificial Complexity Population, Phase I,"
#    my.legend <- "Complexity"
#    # tsps no error bars
#    pdf(file.path(my.sub, "tsps_no_error.pdf"),
#            title=paste(my.title, "TSPs without Errorbars"),
#            width=9, height=6)
#    cat(plot_add_tsps_no_error(my.comp.flow.tsp, my.legend))
#    dev.off()
#    # plot everything for complexity after phase 2
#    my.sub <- file.path(my.dest, "complexity", "final")
#    my.title <- "Artificial Complexity Population, Phase II,"
#    # tsps no error bars
#    pdf(file.path(my.sub, "tsps_no_error.pdf"),
#            title=paste(my.title, "TSPs without Errorbars"))
#    cat(plot_tsps_no_error(my.comp.final.tsp, my.legend))
#    dev.off()
#    # subset of tsps
#    pdf(file.path(my.sub, "tsps_subset.pdf"),
#            title=paste(my.title, "TSPs"))
#    cat(plot_tsps(subset(my.comp.final.tsp, setup %in% c("2.2", "9.3")), my.legend))
#    dev.off()
    # modularity
    Bild1(my.dest)
    Bild2(my.dest)
    Bild3(my.dest)
    Bild4a(my.dest)
    Bild4b(my.dest)
    Bild5(my.dest)
    Bild6(my.dest)
    Bild7(my.dest)
    Bild8(my.dest)
    Bild9(my.dest)
    Bild10(my.dest)
    Bild11(my.dest)
    Bild12(my.dest)
    Bild13(my.dest)
    Bild14(my.dest)
    return("OK")
}

Bild1 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- subset(my.comp.total, type == "Link Robust")
    tmp <- probability_distributions(tmp, c("setup", "phase"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + facet_grid("phase ~ .", scales="free_y")
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S18A_S19A.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase I & II",
            "Link Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild2 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- probability_distributions(my.comp.flow, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "4A.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase I",
            "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild3 <- function(my.dest)
{
    cat(NULL)
    cat("tsps")
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- ddply(my.comp.flow.tsp, c("setup", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE), uncrt=sd(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
            alpha=my.alpha)
    my.plot <- my.plot + geom_errorbar(mapping=aes(ymin=zscore - uncrt,
            ymax=zscore + uncrt, colour=setup), alpha=my.alpha, legend=F)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_brewer(name=my.legend, palette=my.brewer)
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")

    my.file <- file.path(my.dest, "additional", "5A.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase I",
            "TSPs", sep=", ")
    pdf(my.file, title=my.title, width=9, height=5)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild4a <- function(my.dest)
{
    cat(NULL)
    cat("tsps")
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- subset(my.comp.total.tsp, type == "Link Robust")
    tmp <- ddply(tmp, c("phase", "setup", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE), uncrt=sd(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + facet_grid("phase ~ .", scales="free_y")
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
            alpha=my.alpha)
    my.plot <- my.plot + geom_errorbar(mapping=aes(ymin=zscore - uncrt,
            ymax=zscore + uncrt, colour=setup), alpha=my.alpha, legend=F)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_brewer(name=my.legend, palette=my.brewer)
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")

    my.file <- file.path(my.dest, "additional", "5A_6A.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase I & II",
                      "Link Robust", "TSPs", sep=", ")
    pdf(my.file, title=my.title, width=9, height=10)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild4b <- function(my.dest)
{
    cat(NULL)
    cat("tsps")
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final.tsp, type == "Link Robust")
    tmp <- ddply(tmp, c("setup", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE), uncrt=sd(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
            alpha=my.alpha)
    my.plot <- my.plot + geom_errorbar(mapping=aes(ymin=zscore - uncrt,
            ymax=zscore + uncrt, colour=setup), alpha=my.alpha, legend=F)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_brewer(name=my.legend, palette=my.brewer)
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")

    my.file <- file.path(my.dest, "additional", "S13A.pdf")
    my.title <- paste("Parameter Setups", "Phase II", "Link Robust", "TSPs",
                      sep=", ")
    pdf(my.file, title=my.title, width=9, height=5)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild5 <- function(my.dest)
{
    cat(NULL)
    cat("degree correlation")
    my.alpha <- 1
    binw <- 0.02
    my.legend <- "Complexity"

    tmp <- subset(my.comp.total, type == "Node Robust")
    tmp <- probability_distributions(tmp, c("phase", "setup"),
            "degree_correlation", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + facet_grid("phase ~ .", scales="free_y")
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(r)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Degree Correlation ", italic(r))))
    my.plot <- my.plot + coord_cartesian(xlim=c(-1, 1), wise=T)

    my.file <- file.path(my.dest, "additional", "S24B_S25B.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase I & II",
            "Node Robust", "Degree Correlation", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild6 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- subset(my.comp.final, type == "Link Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S19A.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase II",
            "Link Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild7 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- subset(my.comp.final, type == "Node Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S19B.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase II",
            "Node Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild8 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Complexity"

    tmp <- subset(my.comp.final, type == "Noise Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S19C.pdf")
    my.title <- paste("Artificial Complexity Population", "Phase II",
            "Noise Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild9 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Link Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S4A.pdf")
    my.title <- paste("Parameter Setups", "Phase II",
            "Link Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild10 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Node Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S4B.pdf")
    my.title <- paste("Parameter Setups", "Phase II",
            "Node Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild11 <- function(my.dest)
{
    cat(NULL)
    cat("spectral modularity")
    binw <- 0.01
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Noise Robust")
    tmp <- probability_distributions(tmp, c("setup"),
            "spectral_modularity", binw)
    my.plot <- ggplot(tmp, aes(x=x, y=y))
    my.plot <- plot_distribution(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Probability ", P(italic(Q)))))
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Modularity ", italic(Q))))
    my.plot <- my.plot + coord_cartesian(xlim=c(0, 1), wise=TRUE)

    my.file <- file.path(my.dest, "additional", "S4C.pdf")
    my.title <- paste("Parameter Setups", "Phase II",
            "Noise Robust", "Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild12 <- function(my.dest)
{
    cat(NULL)
    cat("complexity vs modularity")
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Link Robust")
    my.plot <- ggplot(tmp, aes(x=scalar_complexity, y=spectral_modularity))
    my.plot <- plot_scatter(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Scalar Complexity ", italic(C))))
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Modularity ", italic(Q))))

    my.file <- file.path(my.dest, "additional", "S2A.pdf")
    my.title <- paste("Parameter Setups", "Phase II", "Link Robust",
            "Scalar Complexity vs Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild13 <- function(my.dest)
{
    cat(NULL)
    cat("complexity vs modularity")
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Node Robust")
    my.plot <- ggplot(tmp, aes(x=scalar_complexity, y=spectral_modularity))
    my.plot <- plot_scatter(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Scalar Complexity ", italic(C))))
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Modularity ", italic(Q))))

    my.file <- file.path(my.dest, "additional", "S2B.pdf")
    my.title <- paste("Parameter Setups", "Phase II", "Node Robust",
            "Scalar Complexity vs Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

Bild14 <- function(my.dest)
{
    cat(NULL)
    cat("complexity vs modularity")
    my.alpha <- 1
    my.legend <- "Parameter Setup"

    tmp <- subset(my.setups.final, type == "Noise Robust")
    my.plot <- ggplot(tmp, aes(x=scalar_complexity, y=spectral_modularity))
    my.plot <- plot_scatter(my.plot, my.legend, my.alpha)
    my.plot <- my.plot + scale_x_continuous(
            name=expression(paste("Scalar Complexity ", italic(C))))
    my.plot <- my.plot + scale_y_continuous(
            name=expression(paste("Modularity ", italic(Q))))

    my.file <- file.path(my.dest, "additional", "S2C.pdf")
    my.title <- paste("Parameter Setups", "Phase II", "Noise Robust",
            "Scalar Complexity vs Spectral (directed) Modularity", sep=", ")
    pdf(my.file, title=my.title)
    cat(my.plot)
    dev.off()
    cat(my.file)
}

plot_add_tsps_no_error <- function(my.df, my.legend, a=1)
{
    cat("tsps")
    cat(paste("alpha =", a))
    tmp <- ddply(my.df, c("setup", "mtf"), summarise,
                 zscore=mean(zscore, na.rm=TRUE))
    my.plot <- ggplot(tmp, aes(x=mtf, y=zscore))
    my.plot <- my.plot + theme_white()
    my.plot <- my.plot + geom_line(mapping=aes(colour=setup, linetype=setup),
            alpha=a)
    my.plot <- my.plot + scale_linetype(name=my.legend)
    my.plot <- my.plot + scale_colour_brewer(name=my.legend, palette=my.brewer)
    my.plot <- my.plot + scale_x_continuous(name="Triad", breaks=1:13)
    my.plot <- my.plot + scale_y_continuous(name="Z-Score")
    my.plot
}


# Old Code ----------------------------------------------------------------


report_correlations <- function(df)
{
    correlate <- factor(c("Sp vs Sc", "Sp vs Bi", "Lo vs Sc", "Lo vs Bi", "Fe vs Sc", "Fe vs Bi"))
    local_df <- data.frame()
    for (typus in levels(df$type)) {
        type_df <- subset(df, type==typus)
        for (name in levels(type_df$setup)) {
            cat(name)
            sub_df <- subset(type_df, setup==name)
            value <- cor(as.numeric(df$spectral_modularity), as.numeric(df$scalar_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Sp vs Sc", type=typus, setup=name))
            value <- cor(as.numeric(df$spectral_modularity), as.numeric(df$binary_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Sp vs Bi", type=typus, setup=name))
            value <- cor(as.numeric(df$louvain_modularity), as.numeric(df$scalar_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Lo vs Sc", type=typus, setup=name))
            value <- cor(as.numeric(df$louvain_modularity), as.numeric(df$binary_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Lo vs Bi", type=typus, setup=name))
            value <- cor(as.numeric(df$mtf_7), as.numeric(df$scalar_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Fe vs Sc", type=typus, setup=name))
            value <- cor(as.numeric(df$mtf_7), as.numeric(df$binary_complexity))
            local_df <- rbind(local_df, data.frame(corr=value, comp="Fe vs Bi", type=typus, setup=name))
        }
    }
    local_df
}

plot_correlations <- function(df)
{
    my.plot <- ggplot(df, aes(x=comp, y=corr, colour=comp))
    my.plot <- my.plot + geom_bar() + stat_identity()
    my.plot <- my.plot + opts(legend.position="none", axis.title.x=theme_blank())
    my.plot <- my.plot + scale_y_reverse(name="Pearson Correlation")
    my.plot
}

