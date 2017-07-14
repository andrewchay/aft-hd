\name{bic.sparsenet}
\alias{bic.sparsenet}
\title{BIC for weighted sparsenet}
\usage{bic.sparsenet(result, X, y)}
\description{Choose the best tuning parameter using BIC criterion.}
\value{A list of results.
\item{betalasso}{The fitted vector of coefficients with LASSO penalty. The length is equal to the number of coefficients in the model.}
\item{betamcp}{The fitted vector of coefficients with MCP penalty. The length is equal to the number of coefficients in the model.}
\item{lasso.lambda}{The best lambda value chosen by BIC criteria when  LASSO is used as the penalty function.}
\item{mcp.lambda}{The best lambda value chosen by BIC criteria when MCP is used as the penalty function.}
\item{mcp.kappa}{The best kappa value chosen by BIC criteria when MCP is used as the penalty function.}
\item{name.lasso}{The names associated with the selected variables by LASSO penalty.}
\item{name.mcp}{The names associated with the selected variables by MCP penalty.}}
\author{Hao Chai <\email{hao.chai@yale.edu}>}
\arguments{\item{result}{a list consists of result from weighted.sparsenet.}
\item{X}{Design matrix of size \eqn{n} by \eqn{p}, where \eqn{n} is the sample size and \eqn{p} is the number of variables. This is the original design matrix.}
\item{y}{Response vector of length \eqn{n}. This is the original response
  variable.}}
\seealso{\code{cv.sparsenet}}
\examples{X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X \%*\% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta, n_kappa = 40, kappa0 = 0.99)
bic.sparsenet(result, X, y)
}

