\name{whiskey}
\alias{whiskey}
\alias{whiskey_brands}
\docType{data}
\title{Survey Data on Brands of Scotch whiskey Consumed}
\description{
  The data set is from Simmons Study of Media and Markets and contains
  the incidence matrix for scotch brands used in last year for those
  households who report consuming scotch.
}
\usage{data("whiskey")}
\format{
  A data frame \code{whiskey} with 484 observations on the following 2
  variables.
  \describe{
    \item{\code{Freq}}{a numeric vector}
    \item{\code{Incidence}}{a matrix with 21 columns}
  }

  Additional information on the brands is contained in the data frame
  \code{whiskey_brands} which is simultaneously loaded. This data frame
  contains 21 observations on the following 3 variables.
  \describe{
    \item{\code{Brand}}{a character vector}
    \item{\code{Type}}{a factor with levels \code{Blend} \code{Single Malt}}
    \item{\code{Bottled}}{a factor with levels \code{Domestic} \code{Foreign}}
  }
}
\details{
  The dataset is taken from the \pkg{bayesm} package. 
}
\source{
  Peter Rossi and Rob McCulloch. (2006). bayesm: Bayesian Inference for
  Marketing/Micro-econometrics. R package version 2.0-8.
  http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html
}

\references{
  Edwards, Y. and G. Allenby (2003), "Multivariate Analysis of
  Multiple Response Data," \emph{JMR} 40, 321-334.
}
\keyword{datasets}
