\name{cv.sparsenet}
\alias{cv.sparsenet}
\title{Cross-validation for weighted sparsenet}
\usage{cv.sparsenet(result, X, y, delta, n_lambda = 100, n_kappa = 1, kappa0 = 1/3, lambda0 = -1, eps = 1e-06, max.iter = 100, method = "include", power = 2, fold = 5, weight = NA, AFT = TRUE, penalty = "MCP", stratified = TRUE)}
\description{Performs k-fold cross validation for MCP penalized regression models over a grid of values for parameter \eqn{\lambda} and \eqn{\kappa}.}
\value{A list of results.
\item{betalasso}{The fitted vector of coefficients with LASSO penalty. The length is equal to the number of coefficients in the model.}
\item{betamcp}{The fitted vector of coefficients with MCP penalty. The length is equal to the number of coefficients in the model.}
\item{betascad}{The fitted vector of coefficients with SCAD penalty. The length is equal to the number of coefficients in the model.}
\item{lasso.lambda}{The best lambda value chosen by BIC criteria when  LASSO is used as the penalty function.}
\item{mcp.lambda}{The best lambda value chosen by BIC criteria when MCP is used as the penalty function.}
\item{mcp.kappa}{The best kappa value chosen by BIC criteria when MCP is used as the penalty function.}
\item{scad.lambda}{The best lambda value chosen by BIC criteria when SCAD is used as the penalty function.}
\item{scad.kappa}{The best kappa value chosen by BIC criteria when SCAD is used as the penalty function.}
\item{name.lasso}{The names associated with the selected variables by LASSO penalty.}
\item{name.mcp}{The names associated with the selected variables by MCP penalty.}}
\author{Hao Chai <\email{hao.chai@yale.edu}>}
\arguments{\item{result}{a list consists of result from weighted.sparsenet.}
\item{X}{\eqn{Design matrix of size \eqn{n} by \eqn{p}, where \eqn{n} is the sample size and \eqn{p} is the number of variables. This is the original design matrix.}}
\item{y}{Response vector of length \eqn{n}. This is the original response variable.}
\item{delta}{Indicator vector of length \eqn{n}. Useful when \code{AFT} is \code{TRUE}. If \eqn{\delta_i} equals 1, \eqn{y_i} is log failure time. If \eqn{\delta_i} equals 0, \eqn{y_i} is log censoring time.}
\item{n_lambda}{number of lambda values in the grid. Default value is 100.}
\item{n_kappa}{number of kappa values in the grid. Default value is 1.}
\item{kappa0}{Maximum value of \eqn{\kappa}. Default value is 1/3.}
\item{lambda0}{Maximum value of \eqn{\lambda}. Default value is -1. If \code{lambda0} = -1, by default the maximum value of \eqn{\lambda} will be the minimum \eqn{\lambda} such that all estimates of \eqn{\beta} are zero.}
\item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than eps. Default is 1e-6.}
\item{max.iter}{Maximum number of iteration. Default is 100.}
\item{method}{Character. Either "exclude" or "include". If \code{method} equals "exclude", after the LASSO step, the function will automatically exclude all the zero estimates from LASSO and pass only the variables associated with non-zero estimates to MCP step. If \code{method} equals "include", then LASSO estimates will serve as a warm start value and all variables will be included in the MCP step.}
\item{power}{Specifies the power of the loss function in cross validation. Default is 2.}
\item{fold}{Number of folds in cross validation. Default is 5.}
\item{weight}{Only useful when \code{AFT} is \code{FALSE}. Specifies the user defined weights associated with observations. Default is \code{NA}.}
\item{AFT}{Logical. If \code{TRUE}, the function uses the accelerated
  failure time model to fit the response variable and \code{weight} will
  be calculated as Kaplan-Meier weights. If \code{FALSE}, \code{weight}
  will be the user specified weights for each observation. Default is
  \code{TRUE}.}
\item{penalty}{The penalty to be added to the objective function. Values could be "MCP" (the default), "SCAD" or "adaptive" for MCP, SCAD and adaptive LASSO penalties, respectively. LASSO result is automatically computed.}
\item{stratified}{Logical. This indicates whether stratified sampling should be used in the cross-validation. If \code{TRUE}, stratified sampling is used. If \code{FALSE}, simple sampling method is used. Default value is \code{TRUE}.}}
\seealso{\code{bic.sparsenet}}
\examples{X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X \%*\% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta, n_kappa = 40, kappa0 = 0.99)
cv.sparsenet(result, X, y, delta, n_kappa = 40, kappa0 = 0.99)
}
