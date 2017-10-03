
\name{bivas}
\alias{bivas}
\title{Fit bivas}
\description{Fit bivas}
\usage{
bivas(y, X, Z = NULL, group, maxIter = 1500, tol = 1e-06, sb2, se2, alpha,
  logodds = NULL, mu, alpha_jk, pi_k, verbose = TRUE, coreNum = 1)
}
\arguments{
\item{y}{response variable. A vector of length p.}

\item{X}{random effect design matrix of n by p.}

\item{Z}{covariate matrix of n by q. Optional. By default, a column of 1's are added.}

\item{group}{a vector of length p, indicating the group of each variable.}

\item{maxIter}{integer. The maximum number of times that EM should be performed. Default is 1500.}

\item{tol}{a number that denotes the tolerance in EM. Default is 1e-8.}

\item{sb2}{a vector of length ns with the elements being the initial value of \sigma_{\beta}^2 for each candidate group prior setting, default value is var(y)/2 for all candidtes.}

\item{se2}{a vector of length ns with the elements being the initial value of \sigma_{e}^2 for each candidate group prior setting, default value is var(y)/2 for all candidates.}

\item{alpha}{a vector of length ns with the elements being the initial value of \alpha, the individual level prior inclusion probability, for each candidate group prior setting, dafault value is 1/log(p+1) for all candidates.}

\item{logodds}{cadidate values of group level prior inclusion probability(pi) in logodds scale, a vector of length ns, default value is seq(-log10(K),log10(K),length.out = 40).}

\item{mu}{initial value of posterior mean of random effects, a p by ns matrix. Default value is 0 for all entries.}

\item{alpha_jk}{initial value of individual level posterior inclusion probility, a p by ns matrix.}

\item{pi_k}{initial value of group level posterior inclusion probility, a K by ns matrix.}

\item{verbose}{Should the values at each iteration to be printed (default=TRUE) or not (FALSE).}

\item{coreNum}{How many physical cores are used in parallel computing, default values is 1 (no parallel).}
}
\details{
  \code{bivas} Fit bivas. It assumes that each column in \code{X} belongs to a uniquely pre-specified group, which is defined in \code{group}. bivas selects both the important groups and the important varibles inside those groups.

  \code{coef.bivas} Returns the coefficient estimates of both random effects and fixec effects.

  \code{predict.bivas} Makes prediction based on fitted \code{bivas} object and new data \code{X}.

  \code{perf.bivas} Evaluate the performance of bivas for a given data using cross-validation.

  \code{getPos.bivas} Returns posterior inclusion probability of both groups and individual variables.

  \code{fdr.bivas} Select important groups and variables by fdr controling.
}
\value{
The function returns a bivas object of S3 class :
  \item{Lq}{
     The values of lower bound of Log-likelihood for each candidate parameter setting.
}
  \item{weight}{
     Normalized weight for each candidate parameter setting.
}
  \item{gamma_pos}{
     Posterior inclusion probability of individual variables.
}
  \item{eta_pos}{
     Posterior inclusion probability of groups.
}
  \item{mu_pos}{
     Posterior mean of individual variables.
}
  \item{cov}{
     Fixed effect estimate for each candidate parameter setting.
}
  \item{sb2}{
     Estimation of parameter \code{sigma2} which is the parameter in spike-slab prior.
}
  \item{se2}{
     Estimation of parameter \code{omega} which measures the proportion of relevant annotations
corresponding to this phenotype.
}
  \item{alpha}{
     Individual level prior probability estimate.
}
  \item{Lq_List}{
     The value of lower bound of Log-likelihood in the main stage.
}
  \item{logodds}{
     The input values of logodds.
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
fit_bivas <- bivas(y=y,X=X,group=g)

}
