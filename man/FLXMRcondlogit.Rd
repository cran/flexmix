\name{FLXMRcondlogit}
\alias{FLXMRcondlogit}
\title{FlexMix Interface to Conditional Logit Models}
\description{
  Model driver for fitting mixtures of conditional logit models.
}
\usage{
FLXMRcondlogit(formula = . ~ ., strata)
}

\arguments{
  \item{formula}{A formula which is interpreted relative to the formula
    specified in the call to \code{\link{flexmix}} using
    \code{\link{update.formula}}. Default is to use the original
    \code{\link{flexmix}} model
    formula.}
  \item{strata}{A formula which is interpreted such that no intercept is
  fitted and which allows to determine the variable indicating which
  observations are from the same stratum.}
}
\details{
  The M-step is performed using \code{coxph.fit}.
}
\value{
  Returns an object of class \code{FLXMRcondlogit} inheriting from \code{FLXMRglm}.
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
\seealso{\code{\link{FLXMRmultinom}}}
\examples{
if (require("Ecdat")) {
  data("Catsup", package = "Ecdat")
  ## To reduce the time needed for the example only a subset is used
  Catsup <- subset(Catsup, id \%in\% 1:100)
  Catsup$experiment <- seq_len(nrow(Catsup))
  vnames <- c("display", "feature", "price")
  Catsup_long <-
    reshape(Catsup,
            idvar = c("id", "experiment"),
            times = c(paste("heinz", c(41, 32, 28), sep = ""),
              "hunts32"),
            timevar = "brand",
            varying = matrix(colnames(Catsup)[2:13], nrow = 3, byrow = TRUE),
            v.names = vnames,
            direction = "long")
  Catsup_long$selected <- with(Catsup_long, choice == brand)
  Catsup_long <- Catsup_long[, c("id", "selected", "experiment", vnames, "brand")]
  Catsup_long$brand <- relevel(factor(Catsup_long$brand), "hunts32")
  set.seed(0808)
  flx1 <- flexmix(selected ~ display + feature + price + brand | id,
                  model = FLXMRcondlogit(strata = ~ experiment), 
                  data = Catsup_long, k = 1)
}
}
\keyword{regression}
\keyword{models}
