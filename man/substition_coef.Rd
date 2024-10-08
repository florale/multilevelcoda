% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{substition_coef}
\alias{substition_coef}
\title{Calculate "coefficients" based on substitutions for each compositional part}
\usage{
substition_coef(object, level = c("between", "within"), h = 10)
}
\arguments{
\item{object}{An object of class \code{brmcoda}.}

\item{level}{A character string specifying the level of the coefficients to be calculated.
Either \dQuote{between} or \dQuote{within}.}

\item{h}{A numeric value specifying the step size for the substitution.}
}
\value{
A data table of results.
}
\description{
Calculate "coefficients" based on substitutions for each compositional part
}
\examples{
\donttest{
if(requireNamespace("cmdstanr")){
sbp2 <- sbp
sbp2[1, ] <- c(-1, 1, -1, -1, -1)
sbp2[2, ] <- c( 1, 0, -1, -1, -1)
sbp2[3, ] <- c( 0, 0,  1, -1, -1)
sbp2[4, ] <- c( 0, 0,  0,  1, -1)

cilr <- complr(data = mcompd, sbp = sbp2, 
  parts = c("TST", "WAKE", "MVPA", "LPA", "SB"),
  idvar = "ID")

m1 <- brmcoda(complr = cilr,
              formula = Stress ~ bilr1 + bilr2 + bilr3 + bilr4 +
                                 wilr1 + wilr2 + wilr3 + wilr4 + (1 | ID),
              chain = 4, iter = 1000, cores = 4L,
              backend = "cmdstanr")
substition_coef(m1, level = "between", h = 10)
substition_coef(m1, level = "within", h = 10)
rm(sbp2, cilr, m1) ## cleanup
}
}
}
