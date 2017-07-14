\name{paths.plot}
\alias{paths.plot}
\title{Plot the solution paths.}
\usage{paths.plot(result, penalty = "MCP", n_kappa = 1)}
\description{This function plots the solution paths of LASSO and MCP solutions.}
\value{Generate the paths plots.}
\author{Hao Chai}
\arguments{\item{result}{a list consists of result from weighted.sparsenet}
\item{penalty}{Character to specify the penalty. Possible values are "MCP", "SCAD", "LASSO" or "adaptive". Specify whether to plot solution paths for MCP, SCAD, LASSO or adaptive LASSO penalties.}
\item{n_kappa}{If \code{penalty} is "MCP", "SCAD" or "adaptive", the function plots the paths corresponding to the \code{n_kappa}th kappa.}}
\examples{X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X \%*\% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta, n_kappa = 40, kappa0 = 0.99)
paths.plot(result, n_kappa = 20)
}
