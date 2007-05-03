%
%  Copyright (C) 2004-2005 Friedrich Leisch
%  $Id: FLXmclust.Rd 1777 2005-08-18 09:19:27Z leisch $
%
\name{plotEll}
\alias{plotEll}
\title{Plot Confidence Ellipses for FLXMCmvnorm Results}
\description{
  Plot confidence ellipses for mixtures of Gaussians fitted using
  \code{\link{FLXMCmvnorm}}.
}
\usage{
plotEll(object, data, which=1:2, project=NULL, points=TRUE,
        eqscale=TRUE, col=NULL, number = TRUE, cex=1.5, numcol="black", 
        pch=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{flexmix} with a fitted
    \code{FLXMCmvnorm} model.}
  \item{data}{The data that were clustered.}
  \item{which}{Index numbers of dimensions of (projected) input space
    to plot.}
  \item{project}{Projection object, currently only the result of
    \code{\link[stats]{prcomp}} is supported.}
  \item{points}{Logical, shall data points be plotted?}
  \item{eqscale}{Logical, plot using \code{\link[MASS]{eqscplot}}?}
  \item{number}{Logical, plot number labels at cluster centers?}
  \item{cex, numcol}{Size and color of number labels.}
  \item{pch, col, \dots}{Graphical parameters.}
}
\author{Friedrich Leisch}
\seealso{\code{\link{FLXMCmvnorm}}}
\keyword{cluster}

