
library(rhdf5)
library(plyr)
library(MASS)


# Data Preparation --------------------------------------------------------


prepare_factors <- function(my.df)
{
    my.df$name <- factor(my.df$name)
    my.df$type <- factor(my.df$type)
    my.df$setup <- factor(my.df$setup)
    return(my.df)
}

reorder_setups <- function(my.df)
{
    my.df$setup <- factor(my.df$setup, levels=c("Standard", "4 Input", "6 Input",
            "10 Input", "12 Input", "2 Activated", "6 Activated", "8 Activated"))
    return(my.df)
}

reorder_complexity <- function(my.df)
{
    tmp <- ddply(my.df, .(setup), colwise(unique, "scalar_complexity"))
    levels(my.df$setup) <- sprintf("%3.1f", round(tmp$scalar_complexity, digits=1))
    my.df$setup <- factor(my.df$setup, levels=c("2.2", "4.0", "9.3", "12.9"))
    return(my.df)
}

load_setups <- function(my.path)
{
    my.h5 <- H5Fopen(my.path)
    # phase 1 setup data
    tmp <- h5read(my.h5, "/setups/phase_1")
    tmp <- prepare_factors(tmp)
    my.setups.p1 <<- reorder_setups(tmp)
    # phase 2 setup data
    tmp <- h5read(my.h5, "/setups/phase_2")
    tmp <- prepare_factors(tmp)
    my.setups.p2 <<- reorder_setups(tmp)
    H5Fclose(my.h5)
    return()
}

load_artificial <- function(my.path)
{
    my.h5 <- H5Fopen(my.path)
    # phase 1 complexity data
    tmp <- h5read(my.h5, "/complexity/phase_1")
    tmp <- prepare_factors(tmp)
    my.comp.p1 <<- reorder_complexity(tmp)
    # phase 2 complexity data
    tmp <- h5read(my.h5, "/complexity/phase_1")
    tmp <- prepare_factors(tmp)
    my.comp.p2 <<- reorder_complexity(tmp)
    H5Fclose(my.h5)
    return()
}

load_all <- function(my.path)
{
    load_setups(my.path)
    load_artificial(my.path)
    return()
}


# Group Prediction --------------------------------------------------------


linear_model <- function(my.df, my.formula)
{
    my.formula <- as.formula(my.formula)
    my.mod <- lm(my.formula, data=my.df)
    my.sum <- summary(my.mod)
    my.res <- as.data.frame(coefficients(my.sum))
#     par(mfrow=c(2, 2))
#     plot(my.mod)
#     par(mfrow=c(1, 1))
    my.res$vars <- rownames(my.res)
    return(my.res)
}

coefficient_agreement <- function(my.ct)
{
    my.N <- sum(my.ct)
    my.f.0 <- sum(diag(my.ct))
    my.rsums <- margin.table(my.ct, 1)
    my.csums <- margin.table(my.ct, 2)
    my.f.c <- sum(my.rsums * my.csums) / my.N
    my.l <- length(my.rsums)
    my.min.marg <- numeric(my.l)
    for (j in 1:my.l) {
        my.min.marg[j] <- min(my.rsums[j], my.csums[j])
    }
    my.f.max <- sum(my.min.marg)
    my.res <- data.frame(p.0=my.f.0 / my.N,
                         p.c=my.f.c / my.N,
                         p.max=my.f.max / my.N,
                         kappa=(my.f.0 - my.f.c) / (my.N - my.f.c),
                         kappa.max=(my.f.max - my.f.c) / (my.N - my.f.c),
                         sigma.kappa=sqrt(my.f.0 * (1 - my.f.0 / my.N)) / (my.N - my.f.c),
                         sigma.kappa.0=sqrt(my.f.c / (my.N * (my.N - my.f.c)))
                         )
    return(my.res)
}

discriminant <- function(my.df, my.formula)
{
    my.formula <- as.formula(my.formula)
    # hackish way of getting the response variable
    my.resp <- all.vars(my.formula)[1]
    my.df[[my.resp]] <- factor(my.df[[my.resp]])
    my.disc <- my.df[[my.resp]]
    for (my.lev in my.disc) {
        if (sum(as.integer(my.disc == my.lev)) < 2) {
            return(data.frame(p.0=NA,
                              p.c=NA,
                              p.max=NA,
                              kappa=NA,
                              kappa.max=NA,
                              sigma.kappa=NA,
                              sigma.kappa.0=NA)
            )
        }
    }
    # lda here could be replaced with another type of analysis
    my.lda <- lda(my.formula, data=my.df, method="mle")
    my.pred <- predict(my.lda)
    my.ct <- table(my.df[[my.resp]], my.pred$class)
    return(coefficient_agreement(my.ct))
}

differentiate_node_noise <- function(my.df, my.frmla="type ~ mean_overlap")
{
    tmp <- subset(my.df, type %in% c("Node Robust", "Noise Robust"))
    tmp$type <- factor(tmp$type)
    return(do_analysis_by_setup(tmp, discriminant, my.frmla))
}

do_analysis_by_setup <- function(my.df, my.func, my.formula)
{
    return(ddply(my.df, .(setup), function(x) {
        my.res <- my.func(x, my.formula)
        my.res$setup <- unique(x$setup)
        return(my.res)
    }))
}
