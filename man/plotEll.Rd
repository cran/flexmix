%
%  Copyright (C) 2004-2015 Friedrich Leisch and Bettina Gruen
%  $Id: plotEll.Rd 5008 2015-01-13 20:30:25Z gruen $
%
\name{plotEll}
\alias{plotEll}
\title{Plot Confidence Ellipses for FLXMCmvnorm Results}
\description{
  Plot 50\% and 95\% confidence ellipses for mixtures of Gaussians fitted using
  \code{\link{FLXMCmvnorm}}.
}
\usage{
plotEll(object, data, which = 1:2, model = 1, project = NULL, points = TRUE,
        eqscale = TRUE, col = NULL, number = TRUE, cex = 1.5, numcol = "black", 
        pch = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{flexmix} with a fitted
    \code{FLXMCmvnorm} model.}
  \item{data}{The response variable in a data frame or as a matrix.}
  \item{which}{Index numbers of dimensions of (projected) input space
    to plot.}
  \item{model}{The model (for a multivariate response) that shall be
    plotted.}
  \item{project}{Projection object, currently only the result of
    \code{\link[stats]{prcomp}} is supported.}
  \item{points}{Logical, shall data points be plotted?}
  \item{eqscale}{Logical, plot using \code{\link[MASS]{eqscplot}}?}
  \item{number}{Logical, plot number labels at cluster centers?}
  \item{cex, numcol}{Size and color of number labels.}
  \item{pch, col, \dots}{Graphical parameters.}
}
\author{Friedrich Leisch and Bettina Gruen}
\seealso{\code{\link{FLXMCmvnorm}}}
\keyword{cluster}

