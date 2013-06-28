
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

inverse_complexity <- function(my.df)
{
    my.df <- transform(my.df, scalar_complexity=1 / scalar_complexity,
                       binary_complexity=1 / binary_complexity)
    return(my.df)
}

reorder_complexity <- function(my.df)
{
    tmp <- ddply(my.df, .(setup), colwise(unique, "scalar_complexity"))
    levels(my.df$setup) <- sprintf("%3.1f", round(tmp$scalar_complexity, digits=1))
    my.df$setup <- factor(my.df$setup, levels=c("2.2", "4.0", "9.3", "12.9"))
    return(my.df)
}

rescale_binary <- function(my.df)
{
    return(ddply(my.df, .(setup), function(x) {
            if (unique(x$setup) == "2 Activated") {
                x$binary_complexity <- x$binary_complexity / 2
            }
            else if (unique(x$setup) == "6 Activated") {
                x$binary_complexity <- x$binary_complexity / 6
            }
            else {
                x$binary_complexity <- x$binary_complexity / 4
            }
            x
        })
    )
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
    # inverse of complexity only necessary for old files
    tmp <- inverse_complexity(tmp)
    my.comp.p1 <<- reorder_complexity(tmp)
    # phase 2 complexity data
    tmp <- h5read(my.h5, "/complexity/phase_1")
    tmp <- prepare_factors(tmp)
    # inverse of complexity only necessary for old files
    tmp <- inverse_complexity(tmp)
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
    my.mod <- lm(my.formula, data=my.df)
#     modularity <- effect("scalar_complexity", my.mod)
#     return(modularity)
    return(my.mod)
}

coefficient_agreement <- function(my.ct)
{
    my.N <- sum(my.ct)
    my.p0 <- sum(diag(my.ct))
    my.rsums <- margin.table(my.ct, 1)
    my.csums <- margin.table(my.ct, 2)
    my.pc <- sum(my.rsums * my.csums) / my.N
    my.min.marg <- numeric(length(my.rsums))
    for (j in 1:length(my.min.marg)) {
        my.min.marg[j] <- min(my.rsums[j], my.csums[j])
    }
    my.p.max <- sum(my.min.marg)
    cat(paste("p0 = ", my.p0, "\n"))
    cat(paste("pc = ", my.pc, "\n"))
    cat(paste("pM = ", my.p.max, "\n"))
    my.kappa <- (my.p0 - my.pc) / (1 - my.pc)
    cat(paste("kappa = ", my.kappa, "\n"))
    my.kappa.max <- (my.p.max - my.pc) / (1 - my.pc)
    cat(paste("kappa max = ", my.kappa.max, "\n"))
    return()
}

linear_discriminant_analysis <- function(my.df, my.frmla="type ~ mean_overlap")
{
    my.frm <- formula(my.frmla)
    my.base <- subset(my.df, type %in% c("Node Robust", "Noise Robust"))
    my.base$type <- factor(my.base$type)
    my.lev <- levels(my.base$setup)
    my.len <- length(my.lev)
    my.p0 <- numeric(my.len)
    my.pc <- numeric(my.len)
    my.k <- numeric(my.len)
    my.k.max <- numeric(my.len)
    my.par <- rep("c", my.len)
#     for (my.set in levels(my.base$setup)) {
    for (i in 1:my.len) {
        my.set <- my.lev[i]
#         cat(paste("\n\n", my.set, "\n", sep=""))
        tmp <- subset(my.base, setup == my.set)
        tmp$setup <- factor(tmp$setup)
#         cat("t.test")
        my1 <- subset(tmp, type == "Node Robust")
        my2 <- subset(tmp, type == "Noise Robust")
#         my.t <- t.test(my1$mean_overlap, my2$mean_overlap)
#         print(my.t)
#         cat("lda\n")
        if (nrow(my1) < 2 | nrow(my2) < 2) {
            my.p0[i] <- NA
            my.pc[i] <- NA
            my.k[i] <- NA
            my.k.max[i] <- NA
            my.par[i] <- my.set
            next
        }
        my.lda <- lda(my.frm, data=tmp, method="mle")
#         print(plot(my.lda))
        my.pred <- predict(my.lda)
#         print(my.lda)
#         print(my.lda$svd)
        my.ct <- table(tmp$type, my.pred$class)
        # elements are normed by total sum of elements
        my.N <- sum(my.ct)
        # sum of diagonal of matrix
        my.f.agree <- sum(diag(my.ct))
        # chance agreement is based on association of marginals
        my.rsums <- margin.table(my.ct, 1)
        my.csums <- margin.table(my.ct, 2)
        my.f.chance <- sum(my.rsums * my.csums) / my.N
        my.min.marg <- numeric(length(my.rsums))
        for (j in 1:length(my.min.marg)) {
            my.min.marg[j] <- min(my.rsums[j], my.csums[j])
        }
        my.f.max <- sum(my.min.marg)
        my.p0[i] <- my.f.agree / my.N
        my.pc[i] <- my.f.chance / my.N
        my.k[i] <- (my.f.agree - my.f.chance) / (my.N - my.f.chance)
        my.k.max[i] <- (my.f.max - my.f.chance) / (my.N - my.f.chance)
        my.par[i] <- my.set
#         my.lda <- lda(my.frm, data=tmp, CV=T, method="mve")
#         ct <- table(tmp$type, my.lda$class)
#         diag(prop.table(ct, 1))
        # total percent correct
#         cat(sum(diag(prop.table(ct))))
#         print("qda")
#         my.qda <- qda(type ~ spectral_modularity + mean_overlap, data=tmp)
#         my.pred <- predict(my.qda)
#         print(table(tmp$type, my.pred$class))
#         my.qda <- qda(type ~ spectral_modularity + mean_overlap, data=tmp, CV=T)
#         ct <- table(tmp$type, my.qda$class)
#         diag(prop.table(ct, 1))
#         # total percent correct
#         print(sum(diag(prop.table(ct))))
    }
    return(data.frame(P0=my.p0, Pc=my.pc, kappa=my.k, kappa_max=my.k.max,
                      setup=my.par))
}
