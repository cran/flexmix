\name{betablocker}
\alias{betablocker}
\docType{data}
\title{Clinical trial of beta-blockers}
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
  G. McLachlan and D. Peel (2000): Finite Mixture Models.
  John Wiley and Sons Inc.
  \url{http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html}
}

\references{
  M. Aitkin (1999): Meta-analysis by random effect modelling in
  generalized linear models.
  Statistics in medicine 18, pages 2343-2351.

  S. Yusuf, R. Peto, J. Lewis, R. Collins and P. Sleight (1985): Beta
  blockade during and after myocardial infarction: an overview of the
  randomized trials.
  Progress in Cardiovascular Diseases 27, pages 335-371.
}
\examples{
data("betablocker")
betaMix <- stepFlexmix(cbind(Deaths, Total-Deaths) ~ 1 | Center,
                       data=betablocker, k=3, nrep=5,
                       model=FLXMRglmfix(family="binomial",
                         fixed=~Treatment)) 
}
\keyword{datasets}
