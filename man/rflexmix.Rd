\name{rflexmix}
\alias{rflexmix}
\alias{rflexmix,flexmix,missing-method}
\alias{rflexmix,FLXdist,numeric-method}
\alias{rflexmix,FLXdist,listOrdata.frame-method}
\alias{rFLXM}
\alias{rFLXM,FLXM,list-method}
\alias{rFLXM,FLXMC,FLXcomponent-method}
\alias{rFLXM,FLXMCmultinom,FLXcomponent-method}
\alias{rFLXM,FLXMCbinom,FLXcomponent-method}
\alias{rFLXM,FLXMRglm,list-method}
\alias{rFLXM,FLXMRglmfix,list-method}
\alias{rFLXM,FLXM,FLXcomponent-method}
\alias{rFLXM,FLXMRglm,FLXcomponent-method}
\title{Random Number Generator for Finite Mixtures}
\description{
  Given a finite mixture model generate random
  numbers from it.
}
\usage{
rflexmix(object, newdata, ...)
}
\arguments{
  \item{object}{A fitted finite mixture model of class \code{flexmix} 
    or an unfitted of class \code{FLXdist}.}
  \item{newdata}{Optionally, a data frame in which to look for variables with
    which to predict or an integer specifying the number of random draws
    for model-based clustering. If omitted, the data to which the model
    was fitted is used.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\details{
  \code{rflexmix} provides the creation of the model matrix for new data
  and the sampling of the cluster memberships. The sampling of the
  component distributions given the classification is done by calling
  \code{rFLXM}. This step has to be provided for the different model
  classes. 
}
\value{
  A list with components
  \item{y}{Random sample}
  \item{group}{Grouping factor}
  \item{class}{Class membership}
}
\author{Bettina Gruen}
\examples{
example(flexmix)
sample <- rflexmix(ex1)
}
\keyword{distribution}


