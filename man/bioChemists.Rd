%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: bioChemists.Rd 4811 2012-05-08 15:35:52Z gruen $
%
\name{bioChemists}
\alias{bioChemists}
\docType{data}
\title{Article production by graduate students in biochemistry Ph.D. programs}
\description{
  A sample of 915 biochemistry graduate students.
}
\usage{
data("bioChemists")
}
\format{
  \describe{
    \item{art}{count of articles produced during last 3 years of Ph.D.}
    \item{fem}{factor indicating gender of student, with levels Men and Women}
    \item{mar}{factor indicating marital status of student, with levels
      Single and Married}
    \item{kid5}{number of children aged 5 or younger}
    \item{phd}{prestige of Ph.D. department}
    \item{ment}{count of articles produced by Ph.D. mentor during last 3 years}
  }
}
\details{
  This data set is taken from package \pkg{pscl} provided by Simon
  Jackman.
}
\source{
  found in Stata format at
  \url{http://www.indiana.edu/~jslsoc/stata/socdata/couart2.dta}
}
\references{
  Long, J. Scott. The origins of sex difference in
  science. \emph{Social Forces}, \bold{68}, 1297--1315, 1990.
  
  Long, J. Scott.  \emph{Regression Models for Categorical and
  Limited Dependent Variables}, 1997. Thousand Oaks, California:
  Sage. 
}
\keyword{datasets}
