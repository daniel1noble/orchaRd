% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_results.R
\name{group_is_valid}
\alias{group_is_valid}
\title{Validate 'group'}
\usage{
group_is_valid(model, group)
}
\arguments{
\item{model}{A meta-analytic model from the \code{metafor} package.}

\item{group}{A character string specifying the name of the grouping variable 
within the model's dataset.}
}
\value{
Logical `TRUE` if the group variable is valid. Otherwise, the function 
throws an error.
}
\description{
Checks if grouping variable is valid within the model's dataset. 
Ensures that the `group` argument is provided, exists as a column 
in the model's data, and is not a numeric continuous variable.
}
\seealso{
\code{\link{mod_results}}
}
\keyword{internal}
