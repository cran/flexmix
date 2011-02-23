%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: patent.Rd 4666 2011-02-23 15:52:35Z gruen $
%
\name{patent}
\alias{patent}
\docType{data}
\title{Patents and R\&D spending}
\description{
  Number of patents, R\&D spending and sales in millions of dollar for 70
  pharmaceutical and biomedical companies in 1976.
}
\usage{data("patent")}
\format{
  A data frame with 70 observations on the following 4 variables.
  \describe{
    \item{Company}{Name of company.}
    \item{Patents}{Number of patents.}
    \item{RDS}{R\&D spending per sales.}
    \item{lgRD}{Logarithmized R\&D spendings (in millions of dollars).}
  }
}
\details{
  The data is taken from the National Bureau of Economic Research R\&D Masterfile.
}
\source{
  P. Wang, I.M. Cockburn and M.L. Puterman (1998): Analysis of Patent
  Data - A Mixed-Poisson-Regression-Model Approach.
  Journal of Business \& Economic Statistics 16 (1), pages 27-41.
}
\references{
  B.H. Hall, C. Cummins, E. Laderman and J. Mundy (1988): The R\&D Master
  File Documentation.
  Technical Working Paper 72, National Bureau of Economic
  Research. Cambridge, MA.
}
\examples{
data("patent")
patentMix <- stepFlexmix(Patents ~ lgRD, k=3,
                      model=FLXMRglm(family="poisson"),
                      concomitant=FLXPmultinom(~RDS),
                      nrep=5, data=patent)
plot(Patents ~ lgRD, data=patent,
     pch=as.character(clusters(patentMix)))
ordering <- order(patent$lgRD)
apply(fitted(patentMix), 2, function(y)
      lines(sort(patent$lgRD), y[ordering]))
}
\keyword{datasets}
