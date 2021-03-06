% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{components}
\alias{components}
\title{Decompose timeseries into trends, seasonal and residuals}
\usage{
components(data, times = NULL)
}
\arguments{
\item{data}{Timeseries to decompose}

\item{times}{Time, expressed in fractional years. If omitted, data are considered monthly}
}
\value{
A list with the decomposed signal
}
\description{
Decompose timeseries into interannual trends, annual and semi-annual signals and residuals
}
\examples{
components( sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01 )
}
