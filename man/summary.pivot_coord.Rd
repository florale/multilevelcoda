% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.pivot_coord}
\alias{summary.pivot_coord}
\title{Create a Summary of a fitted \code{brmsfit} model from a \code{pivot_coord} object}
\usage{
\method{summary}{pivot_coord}(object, digits = 2, ...)
}
\arguments{
\item{object}{An object of class \code{pivot_coord}.}

\item{digits}{A integer value used for number formatting. Default is \code{2}.}

\item{...}{currently ignored.}
}
\value{
A data table of results.
}
\description{
Create a Summary of a fitted \code{brmsfit} model from a \code{pivot_coord} object
}
\examples{
\donttest{
if(requireNamespace("cmdstanr")){
  m <- brmcoda(complr = complr(data = mcompd, sbp = sbp,
                                 parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
                                 idvar = "ID", total = 1440),
  formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
    wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
  chain = 1, iter = 500,
  backend = "cmdstanr")
  
  m_pb <- pivot_coord(m)
  summary(m_pb)
}}
}
