%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: posterior.Rd 3937 2008-03-28 14:56:01Z leisch $
%
\name{relabel}
\alias{relabel}
\alias{relabel,FLXdist,missing-method}
\alias{relabel,FLXdist,character-method}
\alias{relabel,FLXdist,integer-method}
\title{Relabel the Components}

\description{
  The components are sorted by the value of one of the parameters or
  according to an integer vector containing the permutation of the
  numbers from 1 to the number of components.
}

\usage{
relabel(object, by, ...)
\S4method{relabel}{FLXdist,character}(object, by, which = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{"flexmix"}.}
  \item{by}{If a character vector, it needs to be one of \code{"prior"},
    \code{"model"}, \code{"concomitant"} indicating if the parameter
    should be from the component-specific or the concomitant variable
    model. If an integer vector it indicates how the components should
    be sorted.  If missing, the components are sorted by component
    size.}
  \item{which}{Name (or unique substring) of a parameter if \code{by} is
    equal to "model" or "concomitant".}
  \item{\dots}{Currently not used.}
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{methods}
\examples{
    set.seed(123)
    beta <- matrix(1:16, ncol = 4)
    beta
    df1 <- ExLinear(beta, n = 100, sd = .5)
    f1 <- flexmix(y~., data = df1, k = 4)

    ## There was label switching, parameters are not in the same order
    ## as in beta:
    round(parameters(f1))
        
    betas <- rbind(beta, .5)
    betas

    ## This makes no sense:
    summary(abs(as.vector(betas-parameters(f1))))

    ## We relabel the components by sorting the coefficients of x1:
    r1 <- relabel(f1, by = "model", which = "x1")
    round(parameters(r1))

    ## Now we can easily compare the fit with the true parameters:
    summary(abs(as.vector(betas-parameters(r1))))
}
