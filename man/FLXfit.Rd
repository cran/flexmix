\name{FLXfit}
\alias{FLXfit}
\title{Fitter Function for FlexMix Models}
\description{
  This is the basic computing engine called by \code{\link{flexmix}}, it
  should usually not be used directly unless
  by experienced users.
}
\usage{
FLXfit(k = NULL, cluster = NULL, model, control)
}
\arguments{
  \item{k}{Number of clusters (not needed if \code{cluster} is specified).}
  \item{cluster}{Initial cluster assignments of observations at the
    start of the EM algorithm (default is random assignment into \code{k}
    clusters).}
  \item{model}{List of \code{FLXmodel} objects.}
  \item{control}{Object of class \code{FLXcontrol}.}
}
\value{
  Returns an object of class \code{flexmix}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}, \code{\link{flexmix-class}}}
\keyword{regression}
\keyword{cluster}
