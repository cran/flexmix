%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: Mehta.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{Mehta}
\alias{Mehta}
\docType{data}
\title{Mehta Trial}
\description{
  For a 22-centre trial the number of responses and the total number of
  patients is reported for the control group and the group receiving
  a new drug.
}
\usage{data("Mehta")}
\format{
  A data frame with 44 observations on the following 4 variables.
  \describe{
    \item{Response}{Number of responses.}
    \item{Total}{Total number of observations.}
    \item{Drug}{A factor indicating treatment with levels \code{New} and
    \code{Control}.}
    \item{Site}{A factor indicating the site/centre.}
  }
}
\source{
  M. Aitkin. Meta-analysis by random effect modelling in
  generalized linear models.
  \emph{Statistics in Medicine}, \bold{18}, 2343--2351, 1999.
}

\references{
  C.R. Mehta, N.R. Patel and P. Senchaudhuri. Importance sampling
  for estimating exact probabilities in permutational inference.
  \emph{Journal of the American Statistical Association}, \emph{83},
  999--1005, 1988.
}
\examples{
data("Mehta", package = "flexmix")
mehtaMix <- initFlexmix(cbind(Response, Total-Response) ~ 1|Site, 
                        data = Mehta, nrep = 5, k = 3,
                        model = FLXMRglmfix(family = "binomial",
                          fixed = ~ Drug), 
                        control = list(minprior = 0.04))
}
\keyword{datasets}
