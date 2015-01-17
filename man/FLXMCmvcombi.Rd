%
%  Copyright (C) 2009 Friedrich Leisch and Bettina Gruen
%  $Id: FLXMCmvcombi.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{FLXMCmvcombi}
\alias{FLXMCmvcombi}
\title{FlexMix Binary and Gaussian Clustering Driver}
\description{
  This is a model driver for \code{\link{flexmix}} implementing
  model-based clustering of a combination of binary and Gaussian data. 
}
\usage{
FLXMCmvcombi(formula = . ~ .)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
}
\details{
  This model driver can be used to cluster mixed-mode binary and
  Gaussian data. It checks which columns of a matrix contain only zero
  and ones, and does the same as \code{\link{FLXMCmvbinary}} for them. For
  the remaining columns of the data matrix independent Gaussian
  distributions are used (same as \code{\link{FLXMCmvnorm}} with
  \code{diagonal = FALSE}. 
  The same could be obtained by creating a corresponding list of two
  models for the respective columns, but \code{FLXMCmvcombi} does a
  better job in reporting parameters.
}
\value{
  \code{FLXMCmvcombi} returns an object of class \code{FLXMC}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}, \code{\link{FLXMCmvbinary}}, \code{\link{FLXMCmvnorm}}}
\keyword{cluster}
\examples{
## create some artificial data
x1 <- cbind(rnorm(300),
            sample(0:1, 300, replace = TRUE, prob = c(0.25, 0.75)))
x2 <- cbind(rnorm(300, mean = 2, sd = 0.5),
            sample(0:1, 300, replace = TRUE, prob = c(0.75, 0.25)))
x <- rbind(x1, x2)

## fit the model
f1 <- flexmix(x ~ 1, k = 2, model = FLXMCmvcombi())
## should be similar to the original parameters
parameters(f1)
table(clusters(f1), rep(1:2, c(300,300)))

## a column with noise should not hurt too much
x <- cbind(x, rnorm(600))
f2 <- flexmix(x ~ 1, k = 2, model = FLXMCmvcombi())
parameters(f2)
table(clusters(f2), rep(1:2, c(300,300)))

}
