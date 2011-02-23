%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: trypanosome.Rd 4666 2011-02-23 15:52:35Z gruen $
%
\name{trypanosome}
\alias{trypanosome}
\docType{data}
\title{Trypanosome}
\description{
   Trypanosome data from a dosage-response analysis to assess the
   proportion of organisms belonging to different populations. It is
   assumed that organisms belonging to different populations are
   indistinguishable other than in terms of their reaction to the
   stimulus.
}
\usage{data("trypanosome")}
\format{
  A data frame with 426 observations on the following 2 variables.
  \describe{
    \item{\code{Dead}}{A logical vector.}
    \item{\code{Dose}}{A numeric vector.}
  }
}
\details{
  The experimental technique involved inspection under the microscope
  of a representative aliquot of a suspension, all organisms appearing
  within two fields of view being classified either alive or
  dead. Hence the total numbers of organisms present at each dose and
  the number showing the quantal response were both random variables.
}
\source{
  R. Ashford and P.J. Walker (1972): Quantal Response Analysis for a
  Mixture of Populations. Biometrics 28, pages 981-988.
}

\references{
  D.A. Follmann and D. Lambert (1989): Generalizing Logistic
  Regression by Nonparametric Mixing. Journal of the American
  Statistical Association 84(405), pages 195-300.
}

\examples{
data("trypanosome")
trypMix <- stepFlexmix(cbind(Dead, 1-Dead) ~ 1, k = 2, 
	               nrep=5, data = trypanosome,
	               model = FLXMRglmfix(family = "binomial",
                           fixed = ~log(Dose)))
}
\keyword{datasets}
