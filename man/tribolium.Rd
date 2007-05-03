\name{tribolium}
\alias{tribolium}
\docType{data}
\title{Tribolium beetles}
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
  P. Wang and M.L. Puterman (1998): Mixed Logistic Regression Models.
  Journal of Agricultural, Biological, and Environmental Statistics 3
  (2), pages 175-200,
}
\examples{
data("tribolium")
tribMix <- stepFlexmix(cbind(Remaining, Total - Remaining) ~ Species, k = 2, 
                   nrep=5, data = tribolium,
                   model = FLXMRglm(family = "binomial"))
}
\keyword{datasets}
