%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: salmonellaTA98.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{salmonellaTA98}
\alias{salmonellaTA98}
\title{Salmonella Reverse Mutagenicity Assay}
\usage{data("salmonellaTA98")}
\description{
  Data on Ames Salmonella reverse mutagenicity assay.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{Dose levels of quinoline.}
    \item{y}{Numbers of revertant colonies of TA98 Salmonella observed on
      each of three replicate plates tested at each of six dose levels of
      quinoline diameter.}
  }
}

\details{
  This data set is taken from package \pkg{dispmod} provided by Luca
  Scrucca.
}

\source{
  Margolin, B.J., Kaplan, N. and Zeiger, E. Statistical
  analysis of the Ames Salmonella/microsome test,
  \emph{Proc. Natl. Acad. Sci. USA}, \bold{76}, 3779--3783, 1981.
}

\references{
  Breslow, N.E. Extra-Poisson variation in log-linear
  models, \emph{Applied Statistics}, \bold{33}, 38--44, 1984.

  Wang, P., Puterman, M.L., Cockburn, I.M., and Le, N.D. Mixed
  Poisson regression models with covariate dependent rates,
  \emph{Biometrics}, \bold{52}, 381--400, 1996.
}

\examples{
data("salmonellaTA98", package = "flexmix")
salmonMix <- initFlexmix(y ~ 1,
                         data = salmonellaTA98, 
                         model = FLXMRglmfix(family = "poisson", 
                           fixed = ~ x + log(x + 10)),                        
                         k = 2, nrep = 5)
salmonMix.pr <- predict(salmonMix, newdata = salmonellaTA98)
plot(y ~ x, data = salmonellaTA98, 
     pch = as.character(clusters(salmonMix)), 
     ylim = range(c(salmonellaTA98$y, unlist(salmonMix.pr))))
for (i in 1:2) lines(salmonellaTA98$x, salmonMix.pr[[i]], lty = i)
}

\keyword{datasets}
