\name{msharpeBlockSize}
\alias{msharpeBlockSize}
\concept{Testing}
\concept{Sharpe}
\concept{Modified Sharpe}
\title{Optimal block length for bootstrap test of difference of modified Sharpe ratios}
\description{%
  Function which computes the optimal block length for testing of the difference of modified Sharpe ratios}
\usage{%
  msharpeBlockSize(x, y, level = 0.90, na.neg = TRUE, control = list(),
                   b.vec = c(1, 3, 6, 10), alpha = 0.05, M = 199, K = 500, b.av = 5, T.start = 50)
}
\arguments{%
  \item{x}{vector (of lenght \eqn{T}) of returns for the first fund. \code{NA} values are allowed.}
  \item{y}{vector (of lenght \eqn{T}) returns for the second fund. \code{NA} values are allowed.}
  \item{level}{modified Value-at-Risk level. Default: \code{level = 0.90}.}
  \item{na.neg}{a logical value indicating whether \code{NA} values should be returned if a negative modified Value-at-Risk is obtained. 
  Default \code{na.neg = TRUE}.}
  \item{control}{control parameters (see *Details*).}
  \item{b.vec}{vector of block to be tested.}
  \item{alpha}{significance level.}
  \item{M}{bootstrap replications.}
  \item{K}{number of cross-validation.}
  \item{b.av}{average block length in the stationary bootstrap.}
  \item{T.start}{starting point for bootstrap}
}
\details{%
  The argument \code{control} is a list that can supply any of the following components:
  
  \describe{%
    \item{\code{type}}{asymptotic approach (\code{type = 1}) or studentized circular bootstrap approach (\code{type = 2}). Default: \code{type = 1}.}
    \item{\code{ttype}}{test based on ratio (\code{type = 1}) or product (\code{type = 2}). Default: \code{type = 2}.}
    \item{\code{hac}}{heteroscedastic-autocorrelation consistent standard errors. Default: \code{hac = FALSE}.}
    \item{\code{minObs}}{minimum number of concordant observations to compute the ratios. Default: \code{minObs = 10}.}
    \item{\code{nBoot}}{number of boostrap replications for computing the p-value. Default: \code{nBoot = 499}.}
    \item{\code{bBoot}}{block length in the circular bootstrap. Default: \code{bBoot = 1}, i.e. iid bootstrap. \code{bBoot = 0} uses optimal block-length.}
    \item{\code{pBoot}}{symmetric p-value (\code{pBoot = 1}) or asymmetric p-value (\code{pBoot = 2}). Default: \code{pBoot = 1}.}
  }
}
\value{%
   The optimal block length.
}
\note{%
  Please cite the package in publications. Use \code{citation("PeerPerformance")}.
}
\author{%
  David Ardia and Kris Boudt.
}
\references{%
  Ardia, D., Boudt, K. (2015).
  Testing equality of modified Sharpe ratios
  \emph{Finance Research Letters} \bold{13}, pp.97--104.

  Ardia, D., Boudt, K. (2015).
  \emph{The Peer Performance of Hedge Funds}.
  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2000901}

  Favre, L., Galeano, J.A. (2002).
  Mean-modified Value-at-Risk with Hedge Funds.
  \emph{The Journal of Alternative Investments} \bold{5}, pp.21--25.

  Gregoriou, G. N., Gueyie, J.-P. (2003).
  Risk-adjusted performance of funds of hedge funds using a modified Sharpe ratio.
  \emph{The Journal of Wealth Management} \bold{Winter}, pp.77--83.

  Ledoit, O., Wolf, M. (2008). 
  Robust performance hypothesis testing with the Sharpe ratio.
  \emph{Journal of Empirical Finance} \bold{15}, pp.850--859.

  Sharpe, W. F. (1994). 
  The Sharpe ratio. 
  \emph{Journal of Portfolio Management} Fall, pp.49--58.
}
\seealso{%
  \code{\link{msharpe}}, \code{\link{msharpeScreening}} and \code{\link{sharpeTesting}}.
}
\keyword{htest}