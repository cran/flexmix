\name{FLXmclust}
\alias{FLXmclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FlexMix Clustering Demo Driver}
\description{
  This is a demo driver for \code{\link{flexmix}} implementing
  model-based clustering of Gaussian data. 
}
\usage{
FLXmclust(formula = . ~ ., diagonal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Only the left-hand side (response) of
    the formula is used. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{diagonal}{If \code{TRUE}, then the covariance matrix of the
    components is restricted to diagonal matrices.}
}
\details{
  This is meant as a demo for FlexMix driver programming, use package
  mclust for real applications. 
}
\value{
  Returns an object of class \code{FLXmodel}.
}
\author{Friedrich Leisch}
\seealso{\code{\link{flexmix}}}
\keyword{cluster}
\examples{
## The function is currently defined as
FLXmclust

## are packages mvtnorm and ellipse available?
owarn = getOption("warn")
options(warn=-1)
have.mvtnorm=require("mvtnorm", quietly=TRUE)
have.ellipse=require("ellipse", quietly=TRUE)
options(warn=owarn)


if(have.mvtnorm){
  require("MASS")

  ## generate some artificial data
  x <- rbind(rmvnorm(100, mean=rep(0,2)),
             rmvnorm(100, mean=c(8,0), sigma=diag(1:2)),
             rmvnorm(150, mean=c(-2,6), sigma=diag(2:1)),
             rmvnorm(200, mean=c(4,4), sigma=matrix(c(1,.9,.9,1), 2)))

  eqscplot(x)

  ## define a utility function for visualization of results
  plot.ex.mbc = function(object, data){    
    eqscplot(data, col=object@cluster)
    if(have.ellipse){
      for(k in 1:length(object@components)){
          p = parameters(object, k)
          lines(ellipse(p$cov, centre=p$center), col=k)
      } 
    }
  }

  ## This model is wrong (one component has a non-diagonal cov matrix)
  ex.mbc <- flexmix(x~1, k=4, model=FLXmclust())
  print(ex.mbc)
  plot.ex.mbc(ex.mbc, x)

  ## True model, wrong number of components
  ex.mbc <- flexmix(x~1, k=6, model=FLXmclust(diag=FALSE))  
  print(ex.mbc)

  plot.ex.mbc(ex.mbc, x)

  print(parameters(ex.mbc, component=1))
}else{
  cat("Sorry, I need package mvtnorm for these examples!\n\n")
}
}
