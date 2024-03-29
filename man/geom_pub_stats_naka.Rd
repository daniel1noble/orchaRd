% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pub_bias_plot.R
\name{geom_pub_stats_naka}
\alias{geom_pub_stats_naka}
\title{geom_pub_stats_naka}
\usage{
geom_pub_stats_naka(
  data,
  col = "blue",
  plotadj = -0.05,
  textadj = 0.05,
  branch.size = 1.2,
  trunk.size = 3
)
}
\arguments{
\item{data}{The data frame containing the corrected meta-analytic mean and confidence intervals.}

\item{col}{The colour of the mean and confidence intervals.}

\item{plotadj}{The adjustment to the x-axis position of the mean and confidence intervals.}

\item{textadj}{The adjustment to the y-axis position of the mean and confidence intervals for the text displaying the type of correction.}

\item{branch.size}{Size of the confidence intervals.}

\item{trunk.size}{Size of the mean, or central point.}
}
\value{
A list of ggplot2 objects to be added to the orchard plot.
}
\description{
This function adds a corrected meta-analytic mean, sensu Nakagawa et al. 2022, confidence interval and text annotation to an intercept only orchard plot.
}
\author{
Daniel Noble - daniel.noble@anu.edu.au
}
