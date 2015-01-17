%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: tribolium.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{tribolium}
\alias{tribolium}
\docType{data}
\title{Tribolium Beetles}
\description{
  The data investigates whether the adult Tribolium species Castaneum has
  developed an evolutionary advantage to recognize and avoid eggs of
  their own species while foraging.
}
\usage{data("tribolium")}
\format{
  A data frame with 27 observations on the following 4 variables.
  \describe{
    \item{\code{Remaining}}{A numeric vector.}
    \item{\code{Total}}{A numeric vector.}
    \item{\code{Replicate}}{A factor with levels \code{1}, \code{2}, \code{3}.}
    \item{\code{Species}}{A factor with levels \code{Castaneum} \code{Confusum} \code{Madens}.}
  }
}
\details{
  Beetles of the genus Tribolium are cannibalistic in the sense that
  adults eat the eggs of their own species as well as those of closely
  related species. The experiment isolated a number of adult beetles of
  the same species and presented them with a vial of 150 eggs (50 of
  each type), the eggs being thoroughly mixed to ensure uniformity
  throughout the vial.

  The data gives the consumption data for adult Castaneum species. It
  reports the number of Castaneum, Confusum and Madens eggs,
  respectively, that remain uneaten after two day exposure to the adult
  beetles. Replicates 1, 2, and 3 correspond to different occasions on
  which the experiment was conducted.
}
\source{
  P. Wang and M.L. Puterman. Mixed Logistic Regression Models.
  \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{3}
  (2), 175--200, 1998.
}
\examples{
data("tribolium", package = "flexmix")
tribMix <- initFlexmix(cbind(Remaining, Total - Remaining) ~ Species, 
                   k = 2, nrep = 5, data = tribolium,
                   model = FLXMRglm(family = "binomial"))
}
\keyword{datasets}
