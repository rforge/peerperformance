\name{msharpeTesting}
\alias{msharpeTesting}
\concept{Testing}
\concept{Sharpe}
\concept{Modified Sharpe}
\title{Testing the difference of modified Sharpe ratios}
\description{%
  Function which performs the testing of the difference of modified Sharpe ratios}
\usage{%
  msharpeTesting(x, y, level = 0.90, na.neg = TRUE, control = list())
}
\arguments{%
  \item{x}{vector (of lenght \eqn{T}) of returns for the first fund. \code{NA} values are allowed.}
  \item{y}{vector (of lenght \eqn{T}) returns for the second fund. \code{NA} values are allowed.}
  \item{level}{modified Value-at-Risk level. Default: \code{level = 0.90}.}
  \item{na.neg}{a logical value indicating whether \code{NA} values should be returned if a negative modified Value-at-Risk is obtained. 
  Default \code{na.neg = TRUE}.}
  \item{control}{control parameters (see *Details*).}
}
\details{%
  The modified Sharpe ratio (Favre and Galeano 2002) is one industry standard for measuring the absolute risk adjusted performance
  of hedge funds. This function performs the testing of modified Sharpe ratio difference for two funds using 
  a similar approach than Ledoit and Wolf (2002). See also Gregoriou and Gueyie (2003).
  
  For the testing, only the intersection of non-\code{NA} observations for the two funds are used.

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
  A list with the following components:\cr

  \code{n}: number of non-\code{NA} concordant observations.\cr

  \code{msharpe}: vector (of length 2) of unconditional modified Sharpe ratios.\cr

  \code{dmsharpe}: Modified Sharpe ratios difference.\cr

  \code{tstat}: t-stat of modified Sharpe ratios differences.\cr

  \code{pval}: pvalues of test of modified Sharpe ratios differences.
}
\note{%
  Further details on the methdology with an application to the hedge fund industry is given in 
  in Ardia and Boudt (2014). The file \file{ThePeerPerformanceOfHedgeFunds.txt} in the \file{/doc} package's folder
  allows the reprodution of the steps followed in the article. See also the presentation by Kris Boudt 
  at the R/Finance conference 2012 at \url{http://www.rinfinance.com}.

  Some internal functions where adapted from Wolf's R code.

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
\examples{%
  ## Load the data (randomized data of monthly hedge fund returns)
  data("hfdata")
  x = hfdata[,1]
  y = hfdata[,2]	

  ## Run modified Sharpe testing (asymptotic)
  ctr = list(type = 1)
  out = msharpeTesting(x, y, level = 0.95, control = ctr)
  print(out)

  ## Run modified Sharpe testing (asymptotic hac)
  ctr = list(type = 1, hac = TRUE)
  out = msharpeTesting(x, y, level = 0.95, control = ctr)
  print(out)
  
  ## Run modified Sharpe testing (iid bootstrap)
  set.seed(1234)
  ctr = list(type = 2, nBoot = 250)
  out = msharpeTesting(x, y, level = 0.95, control = ctr)
  print(out)

  ## Run modified Sharpe testing (circular bootstrap)
  set.seed(1234)
  ctr = list(type = 2, nBoot = 250, bBoot = 5)
  out = msharpeTesting(x, y, level = 0.95, control = ctr)
  print(out)
}
\keyword{htest}