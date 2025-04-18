% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.is_model_valid}
\alias{.is_model_valid}
\title{Validate Model Argument}
\usage{
.is_model_valid(model)
}
\arguments{
\item{model}{object. A model object from the metafor package (rma.mv, rma, rma.uni, or robust.rma)}
}
\value{
logical. TRUE if the model is valid, otherwise stops with an error message
}
\description{
Checks if a model argument is valid for meta-analysis:
- It is not missing or NULL
- It is a metafor package object (rma.mv, rma, etc.)
}
\author{
Facundo Decunta - fdecunta@agro.uba.ar
}
\keyword{internal}
