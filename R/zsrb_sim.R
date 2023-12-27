#' Calculate the Simple Standard Regression Based Change Score
#'
#' Calculate the Standard Regression Based Change Score as described by
#' Maassen, Bossema, and Brand (https://doi.org/10.1080/13803390802169059).
#' The input data frame must be in wide format.
#'
#' @param idf Input data frame
#' @param t1 Outcome variable scores at time point 1
#' @param t2 Outcome variable scores at time point 2
#' @param group1 Group variable at time point 1
#' @param group2 Group variable at time point 2
#' @param ref Reference group
#' @param covs One or more covariates that you want to take into account
#' @export

zsrb_sim <-
    function(idf, t1, t2, group1, group2, ref, covs = NULL) {

        ## * Select only columns with the Variables of interest
        vois <- c(t1, t2, group1, group2, covs)
        cdf <- idf[, c(vois)]

        ## * Subset the data frame to select only the control group subjects
        # Control subjects should still be in the control subject group on
        # follow-up. E.g., in a study of Alzheimer's disease, a control subject
        # should not have progressed to another group at the second time point.
        ## ** Count the total number of subjects prior to selection
        nrow_a1 <- nrow(cdf)

        ## ** Filter the data
        cdf <- cdf[which(cdf[[group1]] == ref & cdf[[group2]] == ref), ]

        ## ** Count the total number of subjects after selection
        nrow_a2 <- nrow(cdf)

        ## * Select cases with all data available
        # This means that subjects should have data for the outcome measures
        # at both time points, as well as all covariates.

        ## ** Filter the data
        cdf <- cdf[complete.cases(cdf), ]

        ## ** Count the total number of subjects after selection
        nrow_a3 <- nrow(cdf)

        ## * Calculate the estimated beta weight (b)
        # The estimated beta weight (b) was calculated by dividing the standard
        # deviation of the control group at Time 2 (S2) by the standard
        # deviation of the control group at Time 1 (S1), as in b = S2/S1.
        std1 <- sd(cdf[[t1]])
        std2 <- sd(cdf[[t2]])
        b <- std2 / std1

        ## * Calculate the estimated constant (c)
        # The estimated constant (c) was calculated by finding the product of b
        # and the mean of the control group at Time 1 (M1), and subtracting the
        # value from the mean of the control group at Time 2 (M2), as in c = M2
        # – (bM1).
        mt1 <- mean(cdf[[t1]])
        mt2 <- mean(cdf[[t2]])
        c <- mt2 - (b * mt1)

        ## * Calcualte the estimated standard error of the estimate (see)
        # The estimated standard error of the estimate (SEE) was calculated by:
        # - Find the sum of the S1 squared and S2 squared, as in S12+S22.
        # - Find one minus the control group’s test-retest correlation (r12), as
        #   in 1-r12.
        # - Multiply S12+S22 and 1-r12, as in (S12+S22)(1-r12).
        # - Take the square root of the resulting product, as in SEE =
        #   √((S12+S22)(1-r12)).

        ## ** Calculate the sum of S1-squared and S2-squared
        s1s2 <- std1**2 + std2**2

        ## ** Calculate 1 - test-retest correlation
        r12 <- cor(cdf[[t1]], cdf[[t2]], method = "pearson")
        oneminr12 <- 1 - r12

        ## ** Multiply the previous two outcomes
        see_squared <- s1s2 * oneminr12

        ## ** Calculate the see
        see <- sqrt(see_squared)

        ## * Calculate the estimated predicted score at Time 2
        # The estimated predicted Time 2 score (T2’) was calculated by: T2’ =
        # bT1 + c, where T1 is the observed score at Time 1.
        ovar_est <- paste0(t2, "_est")
        idf[[ovar_est]] <- c + (idf[[t1]] * b)

        ## * Estimate the SRB z-score
        # The estimated SRB z-score was calculated by subtracting the T2' score
        # from the observed Time 2 score (T2), and dividing the result by the
        # estimated SEE, as in SRB = (T2 – T2’) / SEE.

        # This z-score can then be compared to the normal distribution, where
        # -1.645 or less traditionally indicates decline, -1.644 to +1.644
        # traditionally indicates no change, and +1.645 or higher traditionally
        # indicates improvement between assessment timepoints.
        ovar_srb <- paste0(t1, "_", t2, "_zsrb")
        idf[[ovar_srb]] <- (idf[[t2]] - idf[[ovar_est]]) / see

        ## * Return the output data frame and some parameters
        ## ** Estimates
        omat <- c(
            "Total number of subjects in the data set", nrow_a1,
            paste("Control subjects that remain controls at", t2), nrow_a2,
            "Control subjects that have complete data", nrow_a3,
            paste("Control subjects: Mean of", t1), mt1,
            paste("Control subjects: Mean of", t2), mt2,
            paste("Control subjects: Standard Deviation of", t1), std1,
            paste("Control subjects: Standard Deviation of", t2), std2,
            "b (std2 / std1)", b,
            "c mt2 - (b * mt1)", c,
            "std1^2 + std2^2", s1s2,
            "1 - (test-retest correlation)", oneminr12,
            "SEE", see
        )
        omat <- as.data.frame(cbind(omat[c(TRUE, FALSE)], omat[c(FALSE, TRUE)]))
        names(omat) <- c("Variable", "Value")

        ## ** Combine output data frame and estimates
        olist <- list(omat, idf)

        ## ** Return output
        return(olist)

    }
