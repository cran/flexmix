%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: plot-methods.Rd 3686 2007-08-14 10:21:46Z gruen $
%
\name{plot-methods}
\docType{methods}
%%%\alias{plot-methods}
\alias{plot,flexmix,missing-method}
\title{Rootogram of Posterior Probabilities}
\description{
  The \code{plot} method for \code{\link{flexmix-class}} objects gives a
  rootogram or histogram of the posterior probabilities.
  }
}
\usage{
\S4method{plot}{flexmix,missing}(x, y, mark=NULL, markcol="red",
  col = NULL, eps=1e-4, root=TRUE, ylim=TRUE, main=NULL, xlab="",
  ylab="", as.table=TRUE, endpoints=c(-0.4, 1.04),...)
}
\arguments{
  \item{x}{an object of class \code{"flexmix"}}
  \item{y}{not used}
  \item{mark}{integer, mark posteriors of this component}
  \item{markcol}{color used for marking components}
  \item{col}{color used for the bars}
  \item{eps}{posteriors smaller than \code{eps} are ignored}
  \item{root}{if \code{TRUE}, a rootogram of the posterior probabilities
    is drawn, otherwise a standard histogram}
  \item{ylim}{a logical value or a numeric vector of length n2. If
    \code{TRUE}, the y axes of all rootograms are aligned
    to have the same limits, if \code{FALSE} each y axis is scaled
    separately. If a numeric vector is specified it is used as usual.}
  \item{main}{main title of the plot}
  \item{xlab}{label of x-axis}
  \item{ylab}{label of y-axis}
  \item{as.table}{logical that controls the order in which panels
    should be plotted: if 'FALSE' (the default), panels are
    drawn left to right, bottom to top (as in a graph); if
    'TRUE', left to right, top to bottom.}
  \item{endpoints}{vector of length 2 indicating the range of x-values that is
    to be covered by the histogram.  This applies only when
    'breaks' is unspecified. In 'do.breaks', this specifies the interval
    that is to be divided up.}
  \item{...}{further graphical parameters for the lattice function
    histogram}
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
  Friedrich Leisch. FlexMix: A general framework for finite mixture
  models and latent class regression in R. Journal of Statistical
  Software, 11(8), 2004. http://www.jstatsoft.org/v11/i08/

  Jeremy Tantrum, Alejandro Murua and Werner Stuetzle. Assessment and
  pruning of hierarchical model based clustering. Proceedings of the
  9th ACM SIGKDD international conference on Knowledge Discovery and
  Data Mining, pages 197-205. ACM Press, New York, NY, USA, 2003.

  Friedrich Leisch. Exploring the structure of mixture model
  components. In Jaromir Antoch, editor, Compstat 2004 - Proceedings in
  Computational Statistics, pages 1405-1412. Physika Verlag, Heidelberg,
  Germany, 2004. ISBN 3-7908-1554-3.
}
\author{Friedrich Leisch}
\keyword{methods}
\keyword{hplot}
