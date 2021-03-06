% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{datediff}
\alias{datediff}
\title{Dates difference}
\usage{
datediff(ymd1, ymd2)
}
\arguments{
\item{ymd1}{Date in YYYYMMDD format (either string or numeric)}

\item{ymd2}{Date in YYYYMMDD format (either string or numeric)}
}
\value{
The days between ymd1 and ymd2
}
\description{
Calculate difference between two dates
}
\examples{
datediff(20120120,20120110)
}
