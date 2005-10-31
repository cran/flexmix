%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: seizure.Rd 1664 2005-06-13 06:11:03Z leisch $
%
\name{seizure}
\alias{seizure}
\docType{data}
\title{Epileptic Seizure Data}
\description{
  Data from a clinical trial where the effect of intravenous
  gamma-globulin on suppression of epileptic seizures is studied.
  There are daily observations for a period of 140 days on one patient, 
  where the first 27 days are a baseline period without treatment, the
  remaining 113 days are the treatment period.
}
\usage{data(seizure)}
\format{
  A data frame with 140 observations on the following 4 variables.
  \describe{
    \item{Seizures}{a numeric vector, daily counts of epileptic seizures}
    \item{Hours}{a numeric vector, hours of daily parental observation}
    \item{Treatment}{a factor with levels \code{No} and \code{Yes}}
    \item{Day}{a numeric vector}
  }
}
\source{
  P. Wang, M. Puterman, I. Cockburn, and N. Le (1996): Mixed poisson
  regression models with covariate dependent rates.
  Biometrics 52, pages 381-400.
}
\references{
  B. Grün and F. Leisch (2004): Bootstrapping finite mixture models.
  In J. Antoch, editor, Compstat 2004 - Proceedings in Computational
  Statistics, pages 1115-1122. Physika Verlag, Heidelberg, Germany,
  ISBN 3-7908-1554-3.
}
\examples{
data(seizure)
plot(Seizures/Hours~Day, col=as.integer(Treatment),
     pch=as.integer(Treatment), data=seizure)
abline(v=27.5, lty=2, col="grey")
legend(140, 9, c("Baseline", "Treatment"),
       pch=1:2, col=1:2, xjust=1, yjust=1)

set.seed(123)

## The model presented in the Wang et al paper: two components for
## "good" and "bad" days, respectively, each a Poisson GLM with hours of
## parental observation as offset

seizMix <- flexmix(Seizures~Treatment*log(Day),
                   data=seizure, k=2,
                   model=FLXglm(family="poisson", offset=log(seizure$Hours)))

summary(seizMix)
summary(refit(seizMix))

matplot(seizure$Day, fitted(seizMix)/seizure$Hours, type="l",
        add=TRUE, col=3:4)
}
\keyword{datasets}
