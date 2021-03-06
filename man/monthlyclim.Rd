% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{monthlyclim}
\alias{monthlyclim}
\title{Compute monthly climatology}
\usage{
monthlyclim(data, mstop = TRUE)
}
\arguments{
\item{data}{Timeseries}

\item{mstop}{Abort if data do not cover years evenly (default to TRUE)}
}
\value{
The 12-month climatology
}
\description{
Compute monthly climatology from inter-annual monthly data.
}
\examples{
monthlyclim(sin( (1:240)*2*pi/12 )+rnorm(240,sd=0.3)+(1:240)*0.01)
}
