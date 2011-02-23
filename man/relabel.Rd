%
%  Copyright (C) 2004-2011 Friedrich Leisch and Bettina Gruen
%  $Id: posterior.Rd 3937 2008-03-28 14:56:01Z leisch $
%
\name{relabel}
\alias{relabel}
\alias{relabel,flexmix,character-method}
\title{Relabel the Components}

\description{
    The components are sorted by the value of one of the
    parameters, and then one EM step is made initialized with the
    posterior probabilities of the previous fit (because of this
    single EM step estimates may change slightly even if EM converged
    previously).  }

\usage{
relabel(object, by, ...)
\S4method{relabel}{flexmix,character}(object, by, data, ...)
}
\arguments{
  \item{object}{An object of class "flexmix".}
  \item{by}{Name (or unique substring) of a parameter.}
  \item{data}{Data frame or list containing the data.  If missing, an
              attempt is made to recover the data from the call, this
              typically works only when used at the prompt and no
              model contains extra variables in the model-specific
              formula .}
  \item{\dots}{Currently not used.}
}
\author{Friedrich Leisch and Bettina Gruen}
\keyword{methods}
\examples{
    set.seed(123)
    beta <- matrix(1:16, ncol=4)
    beta
    df1 <- ExLinear(beta, n=100, sd=.5)
    f1 <- flexmix(y~., data=df1, k=4)

    ## There was label switching, parameters are not in the same order
    ## as in beta:
    round(parameters(f1))
        
    betas <- rbind(beta, .5)
    betas

    ## This makes no sense:
    summary(abs(as.vector(betas-parameters(f1))))

    ## We relabel the components by sorting the coefficients of x1:
    r1 <- relabel(f1, by="x1")
    round(parameters(r1))

    ## Now we can easily compare the fit with the true parameters:
    summary(abs(as.vector(betas-parameters(r1))))
}
