\name{dmft}
\alias{dmft}
\docType{data}
\title{Dental data}
\description{
   Count data from a dental epidemiological study for evaluation of
   various programs for reducing caries collected among school
   children from an urban area of Belo Horizonte (Brazil).
}
\usage{data("dmft")}
\format{
  A data frame with 797 observations on the following 5 variables.
  \describe{
    \item{End}{Number of decayed, missing or filled teeth at the end of the study.}
    \item{Begin}{Number of decayed, missing or filled teeth at the beginning of the study.}
    \item{Gender}{A factor with levels \code{male} and \code{female}.}
    \item{Ethnic}{A factor with levels \code{brown}, \code{white} and \code{black}.}
    \item{Treatment}{A factor with levels \code{control}, \code{educ}, \code{enrich}, 
    \code{rinse}, \code{hygiene} and \code{all}.}
  }
}
  
\details{
   The aim of the caries prevention study was to compare four methods
   to prevent dental caries. Interventions were carried out according
   to the following scheme: 
   \describe{
      \item{control}{Control group}
      \item{educ}{Oral health education}
      \item{enrich}{Enrichment of the school diet with rice bran}
      \item{rinse}{Mouthwash with 0.2\% sodium floride (NaF) solution}
      \item{hygiene}{Oral hygiene}
      \item{all}{All four methods together}
   }
}

\source{
  D. Boehning, E. Dietz, P. Schlattmann, L. Mendonca and U. Kirchner
  (1999): The zero-inflated Poisson model and the decayed, missing and
  filled teeth index in dental epidemiology. Journal of the Royal
  Statistical Society A 162(2), pages 195-209.
}

\examples{
data("dmft")
dmft.flx <- stepFlexmix(End ~ 1, data = dmft, k = 2,
                        model = FLXMRglmfix(family = "poisson", 
                        fixed = ~ Gender + Ethnic + Treatment))
}
\keyword{datasets}
