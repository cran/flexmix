\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,flexmix,missing-method}
\title{Rootogram of Posterior Probabilities}
\description{
  The \code{plot} method for \code{\link{flexmix-class}} objects gives a
  rootogram or histogram of the posterior probabilities.
  }
}
\usage{
\S4method{plot}{flexmix,missing}(x, eps=1e-4, root=TRUE,
  ylim=TRUE, main=NULL, mfrow=NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{"flexmix"}}
  \item{eps}{posteriors smaller than \code{eps} are ignored}
  \item{root}{if \code{TRUE}, a rootogram of the posterior probabilities
    is drawn, otherwise a standard historgram}
  \item{ylim}{A logical value or a numeric vector of length n2. If
    \code{TRUE}, the y axes of all rootograms are aligned
    to have the same limits, if \code{FALSE} each y axis is scaled
    separately. If a numeric vector is specified it is used as usual.}
  \item{main}{main title of the plot}
  \item{mfrow}{layout of the plot}
  \item{...}{further graphical parameters}
}
\details{
  For each mixture component a rootogram or histogram of the posterior
  probabilities of all observations is 
  drawn. Rootograms are very similar to histograms, the only difference
  is that the height of the bars correspond to square roots of counts
  rather than the counts themselves, hence low counts are more
  visible and peaks less emphasized.

  Usually in each component a lot of observations have posteriors
  close to zero, resulting in a high count for the corresponing
  bin in the rootogram which obscures the information in the other
  bins. To avoid this problem, all probabilities with a posterior below
  \code{eps} are ignored.

  A peak at probability one indicates that a mixture component is well
  seperated from the other components, while no peak at one and/or
  significant mass in the middle of the unit interval indicates overlap
  with other components.
}
\references{
  Jeremy Tantrum, Alejandro Murua and Werner Stuetzle. Assessment and
  pruning of hierarchical model based clustering. Proceedings of the
  9th ACM SIGKDD international conference on Knowledge Discovery and
  Data Mining, pages 197-205. ACM Press, New York, NY, USA, 2003.

  Friedrich Leisch. Exploring the structure of mixture model
  components. In Compstat 2004 -- Proceedings in Computational
  Statistics, 2004. Accepted for publication.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
