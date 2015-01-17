%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: FLXbclust.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{FLXMCmvbinary}
\alias{FLXMCmvbinary}
\alias{FLXbclust}
\title{FlexMix Binary Clustering Driver}
\description{
  This is a model driver for \code{\link{flexmix}} implementing
  model-based clustering of binary data. 
}
\usage{
FLXMCmvbinary(formula = . ~ ., truncated = FALSE)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{truncated}{logical, if \code{TRUE} the observations for the
    pattern with only zeros are missing and the truncated likelihood is
    optimized using an EM-algorithm.}
}
\details{
  This model driver can be used to cluster binary data. The only
  parameter is the column-wise mean of the data, which equals the
  probability of observing a 1.
}
\value{
  \code{FLXMCmvbinary} returns an object of class \code{FLXMC}.
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
