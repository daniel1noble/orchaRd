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

