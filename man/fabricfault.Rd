%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: fabricfault.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{fabricfault}
\alias{fabricfault}
\docType{data}
\title{Fabric Faults}
\description{
  Number of faults in rolls of a textile fabric.
}
\usage{data("fabricfault")}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{Length}{Length of role (m).}
    \item{Faults}{Number of faults.}
  }
}
\source{
  G. McLachlan and D. Peel. \emph{Finite Mixture Models}, 2000,
  John Wiley and Sons Inc.
  \url{http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html}
}
\references{
  A. F. Bissell.
  A Negative Binomial Model with Varying Element Sizes
  \emph{Biometrika}, \bold{59}, 435--441, 1972.
  
  M. Aitkin. A general maximum likelihood analysis of
  overdispersion in generalized linear models.
  \emph{Statistics and Computing}, \bold{6}, 251--262, 1996. 
}
\examples{
data("fabricfault", package = "flexmix")
fabricMix <- initFlexmix(Faults ~ 1, data = fabricfault, k = 2,
                         model = FLXMRglmfix(family = "poisson",
                           fixed = ~ log(Length)), 
                         nrep = 5)
}
\keyword{datasets}
