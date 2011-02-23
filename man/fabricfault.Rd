%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: fabricfault.Rd 4666 2011-02-23 15:52:35Z gruen $
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
  G. McLachlan and D. Peel (2000): Finite Mixture Models.
  John Wiley and Sons Inc.
  \url{http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html}
}
\references{
  A. F. Bissell (1972):
  A Negative Binomial Model with Varying Element Sizes
  Biometrika 59, pages 435-441.
  
  M. Aitkin (1996): A general maximum likelihood analysis of
  overdispersion in generalized linear models.
  Statistics and Computing 6, pages 251-262.
}
\examples{
data("fabricfault")
fabricMix <- stepFlexmix(Faults ~ 1, data=fabricfault, k=2,
                         model=FLXMRglmfix(family="poisson",
                           fixed=~ log(Length)), 
                         nrep=5)
}
\keyword{datasets}
