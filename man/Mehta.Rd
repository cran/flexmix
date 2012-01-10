%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: Mehta.Rd 4701 2011-07-20 08:58:44Z gruen $
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
  M. Aitkin (1999): Meta-analysis by random effect modelling in
  generalized linear models.
  Statistics in Medicine 18, pages 2343-2351.
}

\references{
  C.R. Mehta, N.R. Patel and P. Senchaudhuri (1988): Importance sampling
  for estimating exact probabilities in permutational inference.
  Journal of the American Statistical Association 83, pages 999-1005.
}
\examples{
data("Mehta")
mehtaMix <- stepFlexmix(cbind(Response, Total-Response) ~ 1|Site, 
                     data=Mehta, nrep=5, k=3,
                     model=FLXMRglmfix(family="binomial",
                       fixed=~ Drug), 
                     control=list(minprior=0.04))
}
\keyword{datasets}
