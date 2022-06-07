%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: patent.Rd 5220 2022-04-18 18:05:04Z gruen $
%
\name{patent}
\alias{patent}
\docType{data}
\title{Patents and R&D Spending}
\description{
  Number of patents, R&D spending and sales in millions of dollar for 70
  pharmaceutical and biomedical companies in 1976.
}
\usage{data("patent")}
\format{
  A data frame with 70 observations on the following 4 variables.
  \describe{
    \item{Company}{Name of company.}
    \item{Patents}{Number of patents.}
    \item{RDS}{R&D spending per sales.}
    \item{lgRD}{Logarithmized R&D spendings (in millions of dollars).}
  }
}
\details{
  The data is taken from the National Bureau of Economic Research R&D Masterfile.
}
\source{
  P. Wang, I.M. Cockburn and M.L. Puterman. Analysis of Patent
  Data -- A Mixed-Poisson-Regression-Model Approach.
  \emph{Journal of Business & Economic Statistics}, \bold{16}(1),
  27--41, 1998.
}
\references{
  B.H. Hall, C. Cummins, E. Laderman and J. Mundy. The R&D Master
  File Documentation.
  Technical Working Paper 72, National Bureau of Economic
  Research, 1988. Cambridge, MA.
}
\examples{
data("patent", package = "flexmix")
patentMix <- initFlexmix(Patents ~ lgRD, k = 3,
                         model = FLXMRglm(family = "poisson"),
                         concomitant = FLXPmultinom(~RDS),
                         nrep = 5, data = patent)
plot(Patents ~ lgRD, data = patent,
     pch = as.character(clusters(patentMix)))
ordering <- order(patent$lgRD)
apply(fitted(patentMix), 2, function(y)
      lines(sort(patent$lgRD), y[ordering]))
}
\keyword{datasets}
