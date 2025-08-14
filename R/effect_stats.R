



##----------------------------------------------------##
# Internal helper functions
##----------------------------------------------------##

#' @title .calc.skewness
#' @description Computes the skewness of a numeric vector using small sample correction.
#' @param x A numeric vector.
#' @param output Character, either "est" for estimate or "var" for variance.
#' @return A numeric value representing the skewness estimate or variance.
.calc.skewness <- function(x, output = "est") {
  n <- length(x)
  
  if (output == "est") { # skewness estimate
    (sqrt(n * (n - 1)) / (n - 2)) *
      (((1 / n) * sum((x - mean(x)) ^ 3)) /
         (((1 / n) * sum((x - mean(x)) ^ 2)) ^ (3/2)))
    
  } else if (output == "var") { # skewness sampling variance
    (6 * n * (n - 1)) /
      ((n - 2) * (n + 1) * (n + 3))
  }
}

#' @title .calc.kurtosis
#' @description Computes the excess kurtosis of a numeric vector using small sample correction.
#' @param x A numeric vector.
#' @param output Character, either "est" for estimate or "var" for variance.
#' @return A numeric value representing the kurtosis estimate or variance.
.calc.kurtosis <- function(x, output = "est") {
  n <- length(x)
  
  if (output == "est") { # kurtosis estimate
    ((((n + 1) * n * (n - 1)) / ((n - 2) * (n - 3))) *
       (sum((x - mean(x)) ^ 4) / (sum((x - mean(x)) ^ 2) ^ 2))) -(3 * ((n - 1) ^ 2) / ((n - 2) * (n - 3))) # 
  } else if (output == "var") { # kurtosis sampling variance
    (24 * n * ((n - 1) ^ 2)) /
      ((n - 3) * (n - 2) * (n + 3) * (n + 5))
  }
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

