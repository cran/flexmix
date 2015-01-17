%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: betablocker.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{betablocker}
\alias{betablocker}
\docType{data}
\title{Clinical Trial of Beta-Blockers}
\description{
  22-centre clinical trial of beta-blockers for reducing mortality after
  myocardial infarction.
}
\usage{data("betablocker")}
\format{
  A data frame with 44 observations on the following 4 variables.
  \describe{
    \item{Deaths}{Number of deaths.}
    \item{Total}{Total number of patients.}
    \item{Center}{Number of clinical centre.}
    \item{Treatment}{A factor with levels \code{Control} and \code{Treated}.}
  }
}
  
\source{
  G. McLachlan and D. Peel. \emph{Finite Mixture Models}, 2000.
  John Wiley and Sons Inc.
  \url{http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html}
}

\references{
  M. Aitkin. Meta-analysis by random effect modelling in
  generalized linear models.
  \emph{Statistics in Medicine}, \bold{18}, 2343--2351, 1999.

  S. Yusuf, R. Peto, J. Lewis, R. Collins and P. Sleight. Beta
  blockade during and after myocardial infarction: an overview of the
  randomized trials.
  \emph{Progress in Cardiovascular Diseases}, \bold{27}, 335--371,
  1985. 
}
\examples{
data("betablocker", package = "flexmix")
betaMix <- initFlexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center,
                       data = betablocker, k = 3, nrep = 5,
                       model = FLXMRglmfix(family = "binomial",
                         fixed = ~Treatment)) 
}
\keyword{datasets}
