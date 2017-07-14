\name{convexity}
\alias{convexity}
\title{Compute the local convexity of the objective function.}
\usage{convexity(result, eps = sqrt(.Machine$double.eps),
          eigen.min = 0.001, concave.rm = TRUE)}
\description{\code{convexity} computes the minimum eigen value of the Hessian matrix of the penalized weighted least square regression problem.}
\details{For a specific pair of \eqn{\lambda-\kappa}, the \eqn{(i, j)}th entry of the Hessian matrix is computed by \eqn{x'_iWx_j+\frac{\partial^2\sum_{j=1}^p\rho(|\beta_j|\cdot r_{nj};\lambda, \kappa)}{\partial\beta_i\partial\beta_j}}, where is a diagonal matrix with \eqn{i}th diagonal elements \eqn{\omega_{ni}}.}
\value{A list containing \eqn{\lambda}, \eqn{\kappa} etc.
\item{lambda}{vector of length \code{n_lambda}. It contains all the distinct lambda values on the grid.}
\item{kappa}{vector of length \code{n_kappa}.  It contains all the distinct kappa values on the grid.}
\item{index}{A matrix of 4 columns. If \code{concave.rm} is \code{TRUE},
  the number of rows equals the number of \eqn{\lambda-\kappa} pairs on
  the grid at which the solution is convex. If \code{concave.rm} is
  \code{FALSE}, the number of rows equals \code{n_lambda} \eqn{\times} \code{n_kappa}.}}
\author{Hao Chai <\email{hao.chai@yale.edu}>}
\arguments{\item{result}{a list consists of result from weighted.sparsenet.}
\item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than eps. Default is squre root of machine epsilon.}
\item{eigen.min}{The threshold of the minimum eigen value. When the minimum eigen value is smaller than \code{eigen.min}, the minimum eigen value is considered as zero, hence when plug in the optimizer, the minimization problem is locally non-convex. Default is 1e-3.}
\item{concave.rm}{Logical. If \code{TRUE}, when \eqn{\kappa} goes from 0
  to \code{kappa0}, the concave solutions and its associated
  \eqn{\lambda} values will be removed. If \code{FALSE}, when
  \eqn{\kappa} goes from 0 to \code{kappa0}, all \eqn{\lambda} will be
  reserved regardless the convexity of its associated solution. Default
  is \code{TRUE}.}}
\examples{X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X \%*\% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta, n_kappa = 40, kappa0 = 0.99)
a = convexity(result)
}
