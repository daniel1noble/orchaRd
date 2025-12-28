#' @title moment_effects
#' @description Computes the effect estimate of the difference in skewness and kurtosis between two groups along with associated sampling variance for meta-analysis.
#' @param x1 A numeric vector.
#' @param x2 A numeric vector.
#' @param type Character, either "skew" or "kurt" to specify the type of moment effect.
#' @return A data frame with the difference in moment effects and their variances.
#' @examples
#' \dontrun{
#' set.seed(982)
#' library(moments, PearsonDS)
#' 
#' # Just some random comparisons
#' moment_effects(rnorm(100), rnorm(100), type = "skew")
#' moment_effects(rnorm(40), rnorm(40), type = "skew")
#' moment_effects(rnorm(40), rnorm(40), type = "kurt")
#' moment_effects(rnorm(100), rnorm(100), type = "kurt")
#'
#' # Comparisons with known moment differences
#' m1 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 2)
#' m2 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 4)
#'
#' x1 <- PearsonDS::rpearson(1000, moments = m1)
#' x2 <- PearsonDS::rpearson(1000, moments = m2)
#' moment_effects(x1, x2, type = "skew") # ~0
#' moment_effects(x1, x2, type = "kurt") # ~-2
#' 
#' m1 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 3)
#' m2 <- c(mean = 0, variance = 1, skewness = 1, kurtosis = 3)
#' x1 <- PearsonDS::rpearson(1000, moments = m1)
#' x2 <- PearsonDS::rpearson(1000, moments = m2)
#' moment_effects(x1, x2, type = "skew") # ~-1
#' moment_effects(x1, x2, type = "kurt") # ~0
#' }
#' @export
moment_effects <- function(x1, x2, type = c("skew", "kurt")) {
  type <- match.arg(type)

  switch(type,
         skew = {
           return(data.frame(  d_skew = .calc.skewness(x1) - .calc.skewness(x2),
                             d_skew_v = .sv_jack(x1, type = "skew")$var + .sv_jack(x2, type = "skew")$var))
         },
         kurt = {
           return(data.frame(  d_kurt = .calc.kurtosis(x1) - .calc.kurtosis(x2),
                             d_kurt_v = .sv_jack(x1, type = "kurt")$var + .sv_jack(x2, type = "kurt")$var))
         }
	)
}

#' @title Correlation Difference
#' @description Computes the difference in correlation coefficients between two groups along with associated sampling variance for meta-analysis. The correlation and associated sample size or the data frame of the two variables can be provided
#' @param cor1 Numeric, the correlation coefficient for group 1.
#' @param cor2 Numeric, the correlation coefficient for group 2.
#' @param n1 Integer, the sample size for group 1.
#' @param n2 Integer, the sample size for group 2.
#' @param x1 A numeric vector for group 1.
#' @param x2 A numeric vector for group 2.
#' @return A data frame with the difference in correlations and the sampling variances of the difference.
#' @examples
#' \dontrun{
#' set.seed(982)
#' # Example with known correlations
#' cor_diff(0.5, 0.3, 100, 100)
#' cor_diff(0.2, 0.4, 40, 40)
#' }
#' @export
cor_diff <- function(cor1 = NULL, cor2 = NULL, n1 = NULL, n2 = NULL, x1 = NULL, x2 = NULL) {
  df_mode <- !is.null(x1) && !is.null(x2)

  if (df_mode) {
    if (!is.data.frame(x1) || !is.data.frame(x2))
      stop("If x1 and x2 are provided, both must be data frames with exactly two columns.")
    if (ncol(x1) != 2 || ncol(x2) != 2)
      stop("If x1 and x2 are data frames, they must each have exactly two columns.")

    # correlations from columns 1 and 2, pairwise complete
    cor1 <- cor(x1[[1]], x1[[2]], use = "pairwise.complete.obs")
    cor2 <- cor(x2[[1]], x2[[2]], use = "pairwise.complete.obs")

    # effective n = number of complete pairs
    n1 <- sum(stats::complete.cases(x1[[1]], x1[[2]]))
    n2 <- sum(stats::complete.cases(x2[[1]], x2[[2]]))
  } else {
    # numeric path: require r and n for both groups
    if (is.null(cor1) || is.null(cor2) || is.null(n1) || is.null(n2)) {
      stop("Provide either (cor1, cor2, n1, n2) or two data frames x1 and x2 (each with two columns).")
    }
  }

  # sanity checks
  if (n1 <= 3 || n2 <= 3) stop("Sample sizes must be > 3 for Fisher z variance.")
  if (!is.finite(cor1) || !is.finite(cor2)) stop("Correlations must be finite.")

  # Fisher z and variances
  zr1 <- .r.to.zr(cor1)
  zr2 <- .r.to.zr(cor2)
  v1  <- .zr.variance(n1)  # typically 1/(n1 - 3)
  v2  <- .zr.variance(n2)

  data.frame(zr_diff = zr1 - zr2, 
  		 var_zr_diff = v1 + v2)
}

##----------------------------------------------------##
# Internal helper functions
##----------------------------------------------------##

#' @title .calc.skewness
#' @description Computes the skewness of a numeric vector using small sample correction.
#' @param x A numeric vector.
#' @return A numeric value representing the skewness estimate or variance.
.calc.skewness <- function(x) {
	  n <- length(x)
  
   # skewness estimate. Small sample correction.
    (sqrt(n * (n - 1)) / (n - 2)) *
      (((1 / n) * sum((x - mean(x)) ^ 3)) /
         (((1 / n) * sum((x - mean(x)) ^ 2)) ^ (3/2)))
   
}

#' @title .calc.kurtosis
#' @description Computes the excess kurtosis of a numeric vector using small sample correction.
#' @param x A numeric vector.
#' @return A numeric value representing the kurtosis estimate or variance.
.calc.kurtosis <- function(x) {
  n <- length(x)
# Excess kurtosis estimate. Small sample correction.
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
       (sum((x - mean(x)) ^ 4) / (sum((x - mean(x)) ^ 2) ^ 2))) -(3 * ((n - 1) ^ 2) / ((n - 2) * (n - 3))) 
}

#' @title Jackknife Estimation of Skewness and Kurtosis for computing sampling variances
#' @description Computes jackknife estimate and standard error for skewness and kurtosis.
#' @param x A numeric vector.
#' @param bias.correct Logical, whether to apply bias correction. Default is TRUE.
#' @param return.replicates Logical, whether to return jackknife replicates. Default is FALSE.
#' @param type Character, either "skew" or "kurt" to specify the type of statistic.
#' @return A list with estimate (`est`), bias-corrected estimate (`est_bc`), variance, and standard error.

.sv_jack <- function(x, bias.correct = TRUE,
                        return.replicates = FALSE,
						type = c("skew", "kurt")) {
    
	# Match arg
	type <- match.arg(type)

	# Sample size
    n   <- length(x)

	if(type == "skew"){
    	g1  <- function(z) mean((z - mean(z))^3) /
      					   mean((z - mean(z))^2)^(3/2)
	}

	if(type == "kurt"){
    	g1  <- function(z) mean((z - mean(z))^4) /
      					   mean((z - mean(z))^2)^2 - 3
	}

    ## leave-one-out replicates
    theta_i <- vapply(seq_len(n),
                      function(i) g1(x[-i]),
                      numeric(1))
    
      theta   <- g1(x)                 # full-sample estimate
    theta_bar <- mean(theta_i)
    
    ## jack-knife standard error
    se_jack <- sqrt((n - 1) * mean((theta_i - theta_bar)^2))
    
    ## bias correction
    theta_bc <- if (bias.correct) n * theta - (n - 1) * theta_bar else theta
    
    out <- list(est    = theta,    # Estimate
                est_bc = theta_bc, # Bias-corrected estimate
                var    = se_jack^2, # Sampling Variance
                se     = se_jack)   # Standard error

    if (return.replicates) out$jack <- theta_i
    out
}

#' @title Fisher Z-transformation of Correlation Coefficient
#' @description Converts a Pearson correlation coefficient \eqn{r} to Fisher's Zr using the inverse hyperbolic tangent.
#' @param r A numeric vector of Pearson correlation coefficients.
#' @return The Zr-transformed correlation value(s).
.r.to.zr <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

#' @title Variance of Fisher Z-transformed Correlation
#' @description Computes the approximate sampling variance of Fisher's Zr given a sample size.
#' @param n Sample size (must be greater than 3).
#' @return The variance of the Zr-transformed correlation.
.zr.variance <- function(n) {
  1 / (n - 3)
}

#' @title .MSb
#' @description Computes the between group mean-squares difference.
#' @param x1 Mean of group 1.
#' @param x2 Mean of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @return The between group mean-squares difference.
.MSb = function(x1, x2, n1, n2, paired = FALSE){
  n = n1 
  if(paired == TRUE){
    (n/2) * (x1 - x2)^2
  } else {
    ((n1 * n2) / (n1 + n2)) * (x1 - x2)^2
  }
}

#' @title .MSw
#' @description Computes the within group mean-squares difference.
#' @param sd1 Standard deviation of group 1.
#' @param sd2 Standard deviation of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param paired Logical, whether the samples are paired. Default is FALSE.
#' @return The within group mean-squares difference.
.MSw = function(sd1, sd2, n1, n2, paired = FALSE){
  if(paired == TRUE){
    (sd1^2 + sd2^2) / 2
  } else {
    ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  }
}

#' @title .lnM
#' @description Computes the natural logarithm of the ratio of between and within group variances as a measure of the magnitude of difference.
#' @param x1 Mean of group 1.
#' @param x2 Mean of group 2.
#' @param sd1 Standard deviation of group 1.
#' @param sd2 Standard deviation of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @return The natural logarithm of the ratio of between and within group variances.
.lnM <- function(x1, x2, sd1, sd2, n1, n2){
   n0 = (2*n1*n2)/(n1 + n2)
  sw2 = .MSw(sd1, sd2, n1, n2)
  sb2 = (.MSb(x1, x2, n1, n2) - sw2) / n0
  lnM = log(sqrt(sb2)) - log(sqrt(sw2))
  return(lnM)
}


#' @title .v_lnM
#' @description Computes the sampling variance of the natural logarithm of the ratio of between and within group variances.
#' @param x1 Mean of group 1.
#' @param x2 Mean of group 2.
#' @param sd1 Standard deviation of group 1.
#' @param sd2 Standard deviation of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param r Optional correlation between groups for dependent samples. Default is NULL for independent samples.
#' @return The sampling variance of the natural logarithm of the ratio of between and within group variances.

.v_lnM <- function(x1, x2, sd1, sd2, n1, n2, r = NULL){
       n0 = (2*n1*n2) / (n1 + n2)
    theta = x1-x2
        n = n1
  
  if(is.null(r)){
          MSb = .MSb(x1, x2, n1, n2, paired = FALSE)
          MSw = .MSw(sd1, sd2, n1, n2, paired = FALSE)
        delta = MSb - MSw 
         s2_D = sd1^2 / n1 + sd2^2 / n2
        v_lnM = (1 / (4 * delta^2)) * ((n0^2 / 2)^2 * (2 * s2_D^2 + 4 * s2_D * theta^2) + (MSb^2 / MSw^2) * ((2 * MSw^2) / (n1 + n2 -2)))
  } else {
          MSb = .MSb(x1, x2, n1, n2, paired = TRUE)
          MSw = .MSw(sd1, sd2, n1, n2, paired = TRUE)
        delta = MSb - MSw 
         s2_D = sd1^2 + sd2^2 - 2*r*sd1*sd2
        v_lnM = (1 / (4 * delta^2)) * ((n / 2)^2 * (((2 * s2_D^2)/n^2) + ((4 * s2_D * theta^2)/n)) + (MSb^2 / MSw^2) * ((sd1^4 + sd2^4 + 2 * r^2 * sd1^2 * sd2^2) / 2 * (n-1)))
  }
  return(v_lnM)
}


#' @title SAFE bootstrap : independent samples
#' @description Computes a SAFE bootstrap estimate of the effect size and its variance for independent samples.
#' @param x1bar Mean of group 1.
#' @param x2bar Mean of group 2.
#' @param s1 Standard deviation of group 1.
#' @param s2 Standard deviation of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param min_kept Minimum number of accepted bootstrap samples to keep. Default is 2000.
#' @param chunk_init Initial chunk size for bootstrap sampling. Default is 4000.
#' @param chunk_max Maximum chunk size for bootstrap sampling. Default is 2e6.
#' @param max_draws Maximum total number of bootstrap draws. Default is Inf.
#' @param patience_noaccept Number of consecutive chunks with no accepted samples before stopping. Default is 5.
#' @return A list containing the point estimate, variance, number of kept samples, total draws, number of attempts, and status.
#' @export
safe_lnM_indep <- function(x1bar, x2bar, s1, s2, n1, n2,
                           min_kept   = 2000, 
                           chunk_init = 4000,
                           chunk_max  = 2e6,
                           max_draws  = Inf,
                           patience_noaccept = 5) {
  df1 <- n1 - 1L
  df2 <- n2 - 1L
  h   <- (n1 * n2) / (n1 + n2)
  lnM_star <- numeric(0L)
   total <- 0L 
   kept <- 0L
  attempts <- 0L
  zero_streak <- 0L 
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L
    m1 <- rnorm(chunk, mean = x1bar, sd = s1 / sqrt(n1))
    m2 <- rnorm(chunk, mean = x2bar, sd = s2 / sqrt(n2))
    v1 <- s1^2 * rchisq(chunk, df = df1) / df1
    v2 <- s2^2 * rchisq(chunk, df = df2) / df2
    
    MSB <- h * (m1 - m2)^2
    MSW <- (df1 * v1 + df2 * v2) / (df1 + df2)
    good <- which(MSB > MSW)
    n_good <- length(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      vals <- 0.5 * (log((MSB[good] - MSW[good]) / (2 * h)) - log(MSW[good]))
      lnM_star <- c(lnM_star, vals)
      kept <- length(lnM_star)
    } else {
      zero_streak <- zero_streak + 1L
      if (zero_streak >= patience_noaccept) { status <- "no_usable_draws"; break }
    }
    total <- total + chunk
    
    # FIXED ADAPTIVE CHUNKING
    acc <- if (total > 0) kept / total else 0
    if (is_positive(acc)) {
      # Use min() inside next_needed calculation to prevent integer overflow
      next_needed <- min(chunk_max, ceiling(max(0, min_kept - kept) / acc))
    } else {
      next_needed <- chunk * 2L
    }
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  list(point = if (kept >= 2) mean(lnM_star) else NA_real_,
       var   = if (kept >= 2) var(lnM_star) else NA_real_,
       kept  = kept, total = total, attempts = attempts, status = status)
}

#' @title SAFE bootstrap : dependent samples
#' @description Computes a SAFE bootstrap estimate of the effect size and its variance for dependent samples.
#' @param x1bar Mean of group 1.
#' @param x2bar Mean of group 2.
#' @param sd1 Standard deviation of group 1.
#' @param sd2 Standard deviation of group 2.
#' @param n Sample size (assumed equal for both groups).
#' @param r Correlation between the two groups.
#' @param min_kept Minimum number of accepted bootstrap samples to keep. Default is 2000.
#' @param chunk_init Initial chunk size for bootstrap sampling. Default is 4000.
#' @param chunk_max Maximum chunk size for bootstrap sampling. Default is 2e6.
#' @param max_draws Maximum total number of bootstrap draws. Default is Inf.
#' @param patience_noaccept Number of consecutive chunks with no accepted samples before stopping. Default is 5.
#' @return A list containing the point estimate, variance, number of kept samples, total draws, number of attempts, and status.
#' @export
safe_lnM_dep <- function(x1bar, x2bar, sd1, sd2, n, r,
                         min_kept   = 2000,
                         chunk_init = 4000,
                         chunk_max  = 2e6,
                         max_draws  = Inf,
                         patience_noaccept = 5) {
  df <- n - 1L
   h <- n / 2
  Sig <- matrix(c(sd1^2, r*sd1*sd2, r*sd1*sd2, sd2^2), 2, 2)
  lnM_star <- numeric(0L)
  total <- 0L 
  kept <- 0L
  attempts <- 0L
  zero_streak <- 0L 
  chunk <- as.integer(chunk_init)
  status <- "ok"
  
  while (kept < min_kept && total < max_draws) {
    attempts <- attempts + 1L # Added missing 'attempts' incrementer
    Mu <- mvrnorm(n = chunk, mu = c(x1bar, x2bar), Sigma = Sig / n)
    W  <- rWishart(n = chunk, df = df, Sigma = Sig)
    S11 <- W[1,1,] / df 
    S22 <- W[2,2,] / df
    
       MSB <- h * (Mu[,1] - Mu[,2])^2
       MSW <- (S11 + S22) / 2
      good <- which(MSB > MSW)
    n_good <- length(good)
    
    if (n_good > 0L) {
      zero_streak <- 0L
      vals <- 0.5 * (log((MSB[good] - MSW[good]) / n) - log(MSW[good]))
      lnM_star <- c(lnM_star, vals)
      kept <- length(lnM_star)
    } else {
      zero_streak <- zero_streak + 1L
      if (zero_streak >= patience_noaccept) { status <- "no_usable_draws"; break }
    }
    total <- total + chunk
    
    # FIXED ADAPTIVE CHUNKING
    acc <- if (total > 0) kept / total else 0
    if (is_positive(acc)) {
      next_needed <- min(chunk_max, ceiling(max(0, min_kept - kept) / acc))
    } else {
      next_needed <- chunk * 2L
    }
    chunk <- as.integer(max(chunk_init, min(chunk_max, next_needed)))
  }
  
  list(point = if (kept >= 2) mean(lnM_star) else NA_real_,
       var   = if (kept >= 2) var(lnM_star) else NA_real_,
       kept  = kept, total = total, attempts = attempts, status = status)
}
