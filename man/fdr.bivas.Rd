\name{fdr.bivas}
\alias{fdr.bivas}
\title{Select important groups and variables by fdr controling.}
\description{Perform variable selection on both individual variable level and group level based on fdr controling.}
\usage{
fdr.bivas(object,FDRset=0.1,control="local")
}
\arguments{
\item{object}{\code{bivas} object}

\item{FDRset}{The threshold level. Default is 0.1}

\item{control}{Type of control: should be either "local" or global""}

}
\details{
  \code{fdr.bivas} Select important groups and variables by fdr controling.
}
\value{
The function returns a list of :
  \item{FDR}{
     Important variables selected.
}
  \item{groupFDR}{
     Important groups selected.
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
fdr.bivas(object,FDRset=0.1,control="local")

}
