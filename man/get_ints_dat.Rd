% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pub_bias_plot.R
\name{get_ints_dat}
\alias{get_ints_dat}
\title{get_ints_dat}
\usage{
get_ints_dat(model, type = c("bc", "br"))
}
\arguments{
\item{model}{The rma model object containing the corrected meta-analytic mean and confidence intervals.}

\item{type}{The type of correction to extract the corrected meta-analytic mean and confidence intervals from. "br" (i.e., Bias Robust) for Yang et al. 2023, "bc" (i.e., Bias-Corrected) for Nakagawa et al. 2023.}
}
\value{
A list containing the corrected meta-analytic mean and confidence intervals, and a label for the plot.
}
\description{
This function extracts the corrected meta-analytic mean and confidence intervals from a model object.
}
\author{
Daniel Noble - daniel.noble@anu.edu.au
}
