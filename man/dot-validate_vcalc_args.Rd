% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leave_one_out.R
\name{.validate_vcalc_args}
\alias{.validate_vcalc_args}
\title{Validate Variance-Covariance Calculation Arguments}
\usage{
.validate_vcalc_args(model_data, vcalc_args)
}
\arguments{
\item{model_data}{A data frame containing the variables used in the model.}

\item{vcalc_args}{A list of arguments for the metafor::vcalc function.}
}
\value{
The validated vcalc_args list if all checks pass.
}
\description{
Ensures that the arguments provided for variance-covariance calculation are
valid and refer to existing variables in the model data. Performs checks on
the structure and content of the vcalc_args list.
}
\author{
Facundo Decunta - fdecunta@agro.uba.ar
}
\keyword{internal}
