#' Calculate the Simple Standard Regression Based Change Score
#'
#' Calculate the Standard Regression Based Change Score as described by
#' Maassen, Bossema, and Brand (https://doi.org/10.1080/13803390802169059).
#' The input data frame must be in wide format.
#'
#' @param idf Input data frame
#' @param t1 Outcome variable scores at time point 1
#' @param t2 Outcome variable scores at time point 2
#' @param group Group variable. If group membership changes over time, you can
#' add an array with groups, e.g., c("group_t1", "group_t2")
#' @param ref Reference group. If not set, the first factor level or lowest
#' value will be selected
#' @export

zsrb_sim <-
    function(idf, t1, t2, group, ref = NULL) {

        ## * Select only columns with the Variables of interest
        vois <- c(t1, t2, group)
        cdf <- idf[, c(vois)]

        ## * Subset the data frame to select only the control group subjects
        # Control subjects should still be in the control subject group on
        # follow-up. E.g., in a study of Alzheimer's disease, a control subject
        # should not have progressed to another group at the second time point.
        ## ** Count the total number of subjects prior to selection
        nrow_a1 <- nrow(cdf)

        ## ** Filter the data
        ## *** Figure out the reference category if it not set explicitly
        ## If no reference categorie has been set, then select the first
        ## categorie as the reference catergorie
        if (is.null(ref)) {
            if (class(cdf[[group[1]]]) == "numeric")
                ref <- min(cdf[[group[1]]])
        } else if (class(cdf[[group[1]]]) == "factor") {
            ref <- levels(cdf[[group[1]]])[1]
        }

        ## *** Filter the data for reference subjects only
        for (grn in seq_len(length(group))) {
            cdf <- cdf[cdf[[group[grn]]] == ref, ]
        }

        ## *** Remove rownames with NA in the name
        # The previous command can result in rows with 'NA' in the rowname.
        # These rows are compeltely empty. Remove them.
        cdf <- cdf[!grepl("NA", rownames(cdf)), ]

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
        ovar_srb <- paste0("zsrb_", t1, "_", t2)
        idf[[ovar_srb]] <- (idf[[t2]] - idf[[ovar_est]]) / see

        ## * Return the output data frame and some parameters
        ## ** Edit this output table based on if 1 or 2 group vars were entered
        if (grn == 1) {
            tabletext <-
                paste0("No. of `", ref, "` in the sample")
        } else if (grn == 2) {
            tabletext <-
                paste0("No. of `", ref, "` that remain `", ref, "` at ", t2)
        }

        ## ** Estimates
        omat <- c(
            "Total number of observations in the data set", nrow_a1,
            tabletext, nrow_a2,
            paste0("No. of `", ref, "` that have complete data"), nrow_a3,
            paste0("Mean of `", t1, "` (`", ref, "` only; mt1)"), mt1,
            paste0("Mean of `", t2, "` (`", ref, "` only; mt2)"), mt2,
            paste0("St.Dev of `", t1, "` (`", ref, "` only; std1)"), std1,
            paste0("St.Dev of `", t2, "` (`", ref, "` only; std2)"), std2,
            "Estimted beta weight (b): (std2 / std1)", b,
            "Estimated constant (c): mt2 - (b * mt1)", c,
            "Sum of T1 and T2 variance (s1s2): std1^2 + std2^2", s1s2,
            "Test-rest correlation (r12)", r12,
            "1 - 'test-retest correlation' (oneminr12)", oneminr12,
            "SEE: sqrt(s1s2 * oneminr12)", see
        )
        omat <- as.data.frame(cbind(omat[c(TRUE, FALSE)], omat[c(FALSE, TRUE)]))
        names(omat) <- c("Variable", "Value")

        ## ** Combine output data frame and estimates
        olist <- list(omat, idf)

        ## ** Return output
        return(olist)

    }
