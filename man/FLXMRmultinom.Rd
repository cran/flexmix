\name{FLXMRmultinom}
\alias{FLXMRmultinom}
\title{FlexMix Interface to Multiomial Logit Models}
\description{
  Model driver for fitting mixtures of multinomial logit models.
}
\usage{
FLXMRmultinom(formula = . ~ ., ...)
}
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{\dots}{Additional arguments to be passed to \code{nnet.default}.}
}
\details{
  The M-step is performed using \code{nnet.default}.
}
\value{
  Returns an object of class \code{FLXMRmultinom} inheriting from \code{FLXMRglm}.
}
\references{
  Bettina Gruen and Friedrich Leisch. Identifiability of finite mixtures
  of multinomial logit models with varying and fixed
  effects. \emph{Journal of Classification}, \bold{25}, 225--247. 2008.
}
\author{
  Bettina Gruen
}
\section{Warning}{
  To ensure identifiability repeated measurements are
  necessary. Sufficient conditions are given in Gruen and Leisch (2008).
}
\seealso{\code{\link{FLXMRcondlogit}}}
\keyword{regression}
\keyword{models}
