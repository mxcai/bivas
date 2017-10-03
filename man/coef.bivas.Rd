\name{coef.bivas}
\alias{coef.bivas}
\title{Get coefficients of bivas}
\description{Returns the coefficient estimates of both random effects and fixec effects.}
\usage{
coef.bivas(object, weight=T)
}
\arguments{
\item{object}{\code{bivas} object}

\item{weight}{Whether the returned value should be weighted. Default is TRUE.}

}
\details{
  \code{coef.bivas} Returns the coefficient estimates of both random effects and fixec effects.
}
\value{
The function returns a list of :
  \item{beta}{
     Estimate of random effects.
}
  \item{cov}{
     Estimate of fixed effect (including the intercept).
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mingxuan Cai
}
\seealso{
 \code{\link{coef.bivas}}, \code{\link{predict.bivas}}, \code{\link{perf.bivas}}, \code{\link{getPos.bivas}}, \code{\link{fdr.bivas}}
}
\examples{
coef.bivas(object, weight=T)

}
