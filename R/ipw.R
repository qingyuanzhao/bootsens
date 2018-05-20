#' Fast solution to obtain extrema of IPW estimator in missing data
#' problems
#'
#' @param A Indicator of missingness
#' @param X A matrix of covariates
#' @param Y Outcome
#' @param gamma Sensitivity parameter (log odds ratio)
#' @param fitted.prob Fitted propensity score
#' @param estimand Target estimand, either E[Y] ("all") or E[Y|A=0] ("missing")
#'
#' @return Extrema (an interval) of the IPW estimator.
#'
#' @export
#'
get.extrema <- function(A, X, Y, gamma = 0, fitted.prob, estimand = c("all", "missing")) {

    estimand <- match.arg(estimand, c("all", "missing"))
    c <- as.numeric(estimand == "all")

    fitted.logit <- qlogis(fitted.prob)

    eg <- exp(-fitted.logit)
    Y <- Y[A == 1]
    eg <- eg[A == 1]
    eg <- eg[order(-Y)]
    Y <- Y[order(-Y)]

    ## maximization
    num.each.low <- Y * (c + exp(-gamma) * eg)
    num.each.up <- Y * (c + exp(gamma) * eg)
    num <- c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

    den.each.low <- (c + exp(-gamma) * eg)
    den.each.up <- (c + exp(gamma) * eg)
    den <- c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

    maximum <- max(num / den)
    ## print(den[which.max(num/den)] / n)

    ## minimization
    num <- c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
    den <- c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
    minimum <- min(num / den)
    ## print(den[which.min(num/den)] / n)

    c(minimum, maximum)

}

#' Obtain extrema of IPW estimator for missing data
#'
#' @inheritParams get.extrema
#' @param reg.adjust Should regression adjustment (augmented IPW) be used?
#' @param start Starting values for the propensity score model to warm start \code{glm}.
#'
#' @return Extrema (an interval).
#'
#' @import stats
#' @export
#'
extrema.md <- function(A, X, Y, gamma = 0, estimand = c("all", "missing"), reg.adjust = FALSE, start = NULL) {

    estimand <- match.arg(estimand, c("all", "missing"))
    ps.model <- glm(A ~ X, family = "binomial", start = start)
    fitted.prob <- predict(ps.model, type = "response")
    if (reg.adjust) {
        df <- data.frame(Y = Y, X = X)
        or.model <- lm(Y ~ ., df[A == 1, ])
        Y.fitted <- predict(or.model, df)
    } else {
        Y.fitted <- rep(0, length(Y))
    }

    out <- get.extrema(A, X, Y - Y.fitted, gamma, fitted.prob, estimand)
    out <- out + switch(estimand, all = mean(Y.fitted), missing = mean(Y.fitted[A == 0]))

}

#' Obtain extrema of IPW estimator for observational studies
#'
#' @inheritParams extrema.md
#' @param estimand Either "ate" (average treatment effect) or "att" (average treatment effect on the treated)
#'
#' @return Extrema (an interval).
#'
#' @import stats
#' @export
#'
extrema.os <- function(A, X, Y, gamma = 0, estimand = c("ate", "att"), reg.adjust = FALSE, start = NULL) {

    estimand <- match.arg(estimand, c("ate", "att"))

    if (estimand == "att") {
        if (!is.null(start)) {
            mean(Y[A == 1]) - rev(extrema.md(1 - A, X, Y, gamma, "missing", reg.adjust, - start))
        } else {
            mean(Y[A == 1]) - rev(extrema.md(1 - A, X, Y, gamma, "missing", reg.adjust))
        }
    } else { ## estimand == "ate"
        if (!is.null(start)) {
            extrema.md(A, X, Y, gamma, "all", reg.adjust, start) - rev(extrema.md(1 - A, X, Y, gamma, "all", reg.adjust, - start))
        } else {
            extrema.md(A, X, Y, gamma, "all", reg.adjust) - rev(extrema.md(1 - A, X, Y, gamma, "all", reg.adjust))
        }
    }

}


#' Sensitivity analysis for missing data
#'
#' @inheritParams extrema.md
#' @param alpha Significance level
#' @param parallel Should parallel computing be used?
#' @param B Number of Bootstrap resamples.
#' @param warm.start Warm start the refitting of propensity score model (doesn't seem to help).
#'
#' @return A (1 - alpha) confidence interval.
#'
#' @import parallel
#' @export
#'
bootsens.md <- function(A, X, Y, gamma = 0, alpha = 0.05, estimand = c("all", "missing"), reg.adjust = FALSE, parallel = TRUE, B = 1000, warm.start = FALSE) {

    estimand <- match.arg(estimand, c("all", "missing"))

    no.cores <- ifelse(parallel, detectCores(), 1)
    n <- length(A)

    if (warm.start) {
        start <- glm(A ~ X, family = "binomial")$coefs
    } else {
        start <- NULL
    }

    out <- mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE);
        res <- tryCatch(extrema.md(A[s], X[s, ], Y[s], gamma, estimand, reg.adjust, start),
                        error = function(e) {print(e)});
        res},
        mc.cores = no.cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))
}

#' Sensitivity analysis for observational studies
#'
#' @inheritParams bootsens.md
#' @param estimand Either "ate" (average treatment effect) or "att" (average treatment effect on the treated)
#'
#' @return A (1 - alpha) confidence interval.
#'
#' @import stats parallel
#' @export
#'
#' @examples
#'
#' \donttest{
#' ## Fish consumption and blood mercury dataset
#' require(CrossScreening)
#' data(nhanes.fish)
#' A <- as.numeric(nhanes.fish$fish.level == "high")
#' X <- nhanes.fish[, c("gender", "age", "income", "income.missing",
#' "race", "education", "smoking.ever", "smoking.now")]
#' X$race <- factor(X$race)
#' X1 <- model.matrix(~ . - 1, X)
#' Y <- log2(nhanes.fish$o.LBXTHG)
#'
#' ## Assuming no unmeasure confounder (i.e. gamma = 0)
#' extrema.os(A, X1, Y) # point estimate
#' bootsens.os(A, X1, Y) # confidence interval
#'
#' ## Sensitivity analysis (gamma = 1)
#' extrema.os(A, X1, Y, gamma = 1) # point estimate
#' bootsens.os(A, X1, Y, gamma = 1) # confidence interval
#' }
#'
bootsens.os <- function(A, X, Y, gamma = 0, alpha = 0.05, estimand = c("ate", "att"), reg.adjust = FALSE, parallel = TRUE, B = 1000, warm.start = FALSE) {

    estimand <- match.arg(estimand, c("ate", "att"))

    no.cores <- ifelse(parallel, detectCores(), 1)
    n <- length(A)

    if (warm.start) {
        start <- glm(A ~ X, family = "binomial")$coefs
    } else {
        start <- NULL
    }

    out <- mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE);
        res <- tryCatch(extrema.os(A[s], X[s, ], Y[s], gamma, estimand, reg.adjust, start),
                        error = function(e) {print(e)});
        res},
        mc.cores = no.cores)
    out <- do.call(rbind, out)

    c(quantile(out[, 1], alpha / 2, na.rm = TRUE), quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE))

}
