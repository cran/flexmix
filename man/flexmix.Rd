\name{flexmix}
\alias{flexmix}
\alias{flexmix,formula,ANY,ANY,ANY,missing-method}
\alias{flexmix,formula,ANY,ANY,ANY,list-method}
\alias{flexmix,formula,ANY,ANY,ANY,FLXmodel-method}
\alias{parameters}
\alias{parameters,flexmix-method}
\alias{predict,flexmix-method}
\alias{show,flexmix-method}
\title{Flexible Mixture Modeling}
\description{
  FlexMix implements a general framework for finite
  mixtures of regression models. Parameter estimation is performed using
  the EM algorithm: the E-step is implemented by \code{flexmix}, while
  the user can specify the M-step.
}
\usage{
flexmix(formula, data = list(), k = NULL, cluster = NULL,
        model=NULL, control = NULL)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit. The details of
    model specification are given below.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{k}{Number of clusters (not needed if \code{cluster} is specified).}
  \item{cluster}{Factor or integer vector with the initial cluster
    assignments of observations at the 
    start of the EM algorithm (default is random assignment into \code{k}
    clusters).}
  \item{model}{Object of \code{FLXmodel} of list of \code{FLXmodel}
    objects. Default is the object returned by calling
    \code{\link{FLXglm}()}.}
  \item{control}{Object of class \code{FLXcontrol} or a named list.}
}
\details{
  FlexMix models are described by objects of class \code{FLXmodel},
  which in turn are created by driver functions like
  \code{\link{FLXglm}} or \code{\link{FLXmclust}}. Multivariate
  responses with independent components can be specified using a
  list of \code{FLXmodel} objects. 
}
\value{
  Returns an object of class \code{flexmix}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix-class}}, \code{\link{FLXglm}},
  \code{\link{FLXmclust}}}
\keyword{regression}
\keyword{cluster}
\examples{

## generate some artificial data with a Gaussian response y and a
## Poisson response yp
x <- runif(200, 0, 10)
mydf <- data.frame(x=x,
                 y=c(5*x[1:100], 40-(x[101:200]-5)^2)+3*rnorm(200),
                 yp=rpois(200, exp(c(2-0.2*x[1:100], 1+0.1*x[101:200]))),
                 class=rep(1:2, c(100,100)))

plot(y~x, data=mydf, col=class)

## mixture of two linear regression models. Note that control parameters
## can be specified as named list and abbreviated if unique.
ex.qpol <- flexmix(y~x+I(x^2), data=mydf, k=2,
                   control=list(verb=5, iter=100))

ex.qpol
ex.qpol@components

## now use both responses simultanously
## first we have a look at the data
opar <- par("mfrow")
par(mfrow=1:2)
plot(y~x, data=mydf, col=class)
plot(yp~x, data=mydf, col=class)
par(mfrow=opar)

## now we fit a model with one Gaussian response and one Poisson
## response. Note that the formulas inside the call to FLXglm are
## relative to the overall model formula.
ex.qpol1 <- flexmix(y~x, data=mydf, k=2,
                    model=list(FLXglm(y~.+I(x^2)), 
                               FLXglm(yp~., family="poisson")))  

ex.qpol1
table(ex.qpol1@cluster, mydf$class)

## for Gaussian responses we get coefficients and standard deviation
parameters(ex.qpol1, component=1, model=1)

## for Poisson response we get only coefficients
parameters(ex.qpol1, component=1, model=2)
}

