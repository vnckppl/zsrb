#' Calculate the Simple or Complex Standard Regression Based Change Score
#'
#' Calculate the Standard Regression Based Change Score as described by
#' Hammers, Suhrie, Dixon, Porter and Duff (2016).
#' https://doi.org/10.1093/arclin/acz076
#' The input data frame must be in wide format.
#'
#' @param idf Input data frame
#' @param t1 Outcome variable scores at time point 1
#' @param t2 Outcome variable scores at time point 2
#' @param group Group variable. If group membership changes over time, you can
#' add an array with groups, e.g., c("group_t1", "group_t2")
#' @param ref Reference group. If not set, the first factor level or lowest
#' value will be selected
#' @param covs One or more covariates that you want to take into account
#' @export

zsrb_com <- function(idf, t1, t2, group, ref = NULL, covs = NULL) {

    ## * Select only columns with the Variables of interest
    vois <- c(t1, t2, group, covs)
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
    # The previous command can result in rows with 'NA' in the rowname. These
    # rows are compeltely empty. Remove them.
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

    ## * Build and run the regression model
    # Build the regression model that estimates the T2 score from the T1
    # score, group variable and the covariates. Note that 'group' is not
    # included here because we have subsetted our data to only include the
    # reference group.
    variables <- paste(c(t1, covs), collapse = " + ")
    my_formula <- paste(t2, "~", variables)
    model <- lm(formula = my_formula, data = cdf)

    ## ** Store the coefficients for easy access
    coeffs <- summary(model)$coefficient

    ## * Itteratively calculate t2' score for all subjects
    # t2' is the estimated score using the intercept and beta values that
    # we extract from the regression model. This process is itterative,
    # because we loop over all covariates. First, let's generate the t2'
    # variable name.
    t2_prime <- paste0(t2, "_prime")

    ## ** Use intercept and T1 score to start the estimate of t2_prime score
    idf[[t2_prime]] <-
        coeffs["(Intercept)", "Estimate"] +
        idf[[t1]] * coeffs[t1, "Estimate"]

    ## ** Loop over covariates and add their effects to t2_prime
    # Only do this if there are covariates
    if (!is.null(covs)) {
        for (covar in covs) {
            ## *** For numeric covariates
            # There is one beta, so simply add the effect
            if (
                class(idf[[covar]]) == "numeric" ||
                class(idf[[covar]]) == "integer")
            {
                idf[[t2_prime]] <- idf[[t2_prime]] +
                    idf[[covar]] * coeffs[covar, "Estimate"]

            ## *** For difftime variables
            # This works the same as numeric variables
            } else if (class(idf[[covar]]) == "difftime") {
                idf[[t2_prime]] <- idf[[t2_prime]] +
                    as.numeric(idf[[covar]]) * coeffs[covar, "Estimate"]

            ## *** For categorical covariates
            # There are one or more betas in this case, so we need to take
            # all of them into account
            } else if (class(idf[[covar]]) == "factor") {
                ## **** Loop over levels of the factor variable
                # We skip the first factor level, because that is the
                # reference level, so it does not have a beta value
                for (fact in levels(idf[[covar]])[2:length(levels(idf[[covar]]))]) {

                    ## ***** Construct the name of this var in the coefs
                    # Construct the name of this variable in the coefficient
                    # output from the linear model. It is just the
                    # combination of the variable name and the factor level.
                    facname <- paste0(covar, fact)

                    ## ***** Add the effect of the covariate to the data
                    # But only do this for the rows where the covariate is
                    # that covariate level (e.g., if there is a beta for
                    # males, then only apply that to the rows with males)
                    idf[[t2_prime]][idf[[covar]] == fact] <-
                        idf[[t2_prime]][idf[[covar]] == fact] +
                        coeffs[facname, "Estimate"]

                }
            }
        }
    }

    ## * Calculate the difference between T2' and T2
    t2_diff <- idf[[t2]] - idf[[t2_prime]]

    ## * Standardize this value by the model's 'Residual Standard Error'
    # The residual standard error in R's lm is what Hammers and Duff refer
    # to as the 'standard error of the estimate of the regression model', or
    # SEE as they derived it from SPSS output.
    see <- sigma(model)

    ## * Now calculate the zSRB
    zsrb_var <- paste("zsrb", t1, t2, sep = "_")
    zsrb_score <- t2_diff / see
    idf[[zsrb_var]] <- zsrb_score

    ## * Return the output data frame and some parameters
    ## ** Estimates
    omat <- c(
        "Total number of subjects in the data set", nrow_a1,
        paste(ref, "subjects that remain", ref, "at", t2), nrow_a2,
        paste(ref, "subjects that have complete data"), nrow_a3,
        "Regression formula", my_formula,
        "Reference group", ref
    )
    omat <- as.data.frame(cbind(omat[c(TRUE, FALSE)], omat[c(FALSE, TRUE)]))
    names(omat) <- c("Variable", "Value")

    ## ** Combine output data frame and estimates
    olist <- list(omat, model, idf)

    ## ** Return output
    return(olist)

}
