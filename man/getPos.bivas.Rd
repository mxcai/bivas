\name{getPos.bivas}
\alias{getPos.bivas}
\title{Evaluate the posterior inclusion probability.}
\description{Calculate the posterior inclusion probability from \code{bivas} object.}
\usage{
getPos.bivas(object, weight=T)
}
\arguments{
\item{object}{\code{bivas} object}

\item{weight}{Whether the probability should be weighted.Default is TRUE.}


}
\details{
  \code{getPos.bivas} Returns posterior inclusion probability of both groups and individual variables.
}
\value{
The function returns a list of :
  \item{var_pos}{
     The posterior inclusion probility of individual variables.
}
  \item{group_pos}{
     The posterior inclusion probility of groups.
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
getPos.bivas(object, weight=T)

}
