%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: seizure.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{seizure}
\alias{seizure}
\docType{data}
\title{Epileptic Seizure Data}
\description{
  Data from a clinical trial where the effect of intravenous
  gamma-globulin on suppression of epileptic seizures is studied.
  Daily observations for a period of 140 days on one patient are given, 
  where the first 27 days are a baseline period without treatment, the
  remaining 113 days are the treatment period.
}
\usage{data("seizure")}
\format{
  A data frame with 140 observations on the following 4 variables.
  \describe{
    \item{Seizures}{A numeric vector, daily counts of epileptic seizures.}
    \item{Hours}{A numeric vector, hours of daily parental observation.}
    \item{Treatment}{A factor with levels \code{No} and \code{Yes}.}
    \item{Day}{A numeric vector.}
  }
}
\source{
  P. Wang, M. Puterman, I. Cockburn, and N. Le. Mixed poisson
  regression models with covariate dependent rates.
  \emph{Biometrics}, \bold{52}, 381--400, 1996.
}
\references{
  B. Gruen and F. Leisch. Bootstrapping finite mixture models.
  In J. Antoch, editor, Compstat 2004--Proceedings in Computational
  Statistics, 1115--1122. Physika Verlag, Heidelberg, Germany, 2004.
  ISBN 3-7908-1554-3.
}
\examples{
data("seizure", package = "flexmix")
plot(Seizures/Hours ~ Day, col = as.integer(Treatment),
     pch = as.integer(Treatment), data = seizure)
abline(v = 27.5, lty = 2, col = "grey")
legend(140, 9, c("Baseline", "Treatment"),
       pch = 1:2, col = 1:2, xjust = 1, yjust = 1)

set.seed(123)

## The model presented in the Wang et al paper: two components for
## "good" and "bad" days, respectively, each a Poisson GLM with hours of
## parental observation as offset

seizMix <- flexmix(Seizures ~ Treatment * log(Day),
                   data = seizure, k = 2,
                   model = FLXMRglm(family = "poisson",
                     offset = log(seizure$Hours)))

summary(seizMix)
summary(refit(seizMix))

matplot(seizure$Day, fitted(seizMix)/seizure$Hours, type = "l",
        add = TRUE, col = 3:4)
}
\keyword{datasets}
