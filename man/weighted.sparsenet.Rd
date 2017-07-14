\name{weighted.sparsenet}
\alias{weighted.sparsenet}
\title{Fit penalized weighted least square regression}
\usage{weighted.sparsenet(X, y, delta = 1, n_lambda = 100, n_kappa = 1, kappa0 = 1/3, lambda0 = -1, eps = 1e-06, max.iter = 1000, method = "include", name.list = NA, weight = NA, AFT = TRUE, penalty = "MCP")}
\description{\code{weighted.sparsenet} is a function to fit the penalized weighted least square regression using coordinate descent algorithm.}
\details{The function \code{weighted.sparsenet} computes the solutions
  to the penalized least square problem. It minimizes \deqn{L_{\lambda,
      \gamma}(\beta_0, \beta) =
    1/2\sum_{i=1}^n\omega_{ni}(Y_i-\beta_0
    -X'_i\beta)^2+\sum\limits_{j=1}^p\rho(r_{nj}|\beta_j|;\lambda,
    \gamma),} where \eqn{\omega_{n1}=\delta_{(1)}/n}, \eqn{\omega_{ni}=\delta_{(i)}/(n-i+1)\prod_{j=1}^{i-1}((n-j)/(n-j+1))^{\delta_{(j)}}} for \eqn{i=2, 3, \ldots, n}, and \eqn{r_{nj}=\sqrt{\sum\limits_{i=1}^n\omega_{ni}(X_{ij}-\bar X_{\omega j})^2}}.
For \eqn{t>0}, when the penalty is LASSO,
\eqn{\rho(t;\lambda)=|t|}. When the penalty is MCP,
\eqn{\rho(t;\lambda,\kappa)=\int_0^t\lambda(1-\frac{\kappa
    x}{\lambda})_+dx}. In the first step, it computes the solutions with LASSO penalty along \eqn{\lambda}. In the second step, it uses the LASSO solution  as a warm start value to compute the solution with MCP penalty. At each point on the \eqn{\lambda - \kappa} surface, the function will return the corresponding solutions. To form the grid, a sequence of values of length \code{n_lambda} is computed, equally spaced on the log scale and is assigned to \eqn{\lambda}, and a sequence of values of length \code{n_kappa} is computed, equally spaced and is assigned to \eqn{\lambda}. The minimum of \eqn{\lambda} is 0 and the maximum of \eqn{\lambda} is \code{lambda0}. The minimum of \eqn{\kappa} is 0 and the maximum of \eqn{\kappa} is \code{kappa0}.}
\value{A list of results.
\item{omega}{vector of length \eqn{n}. It contains the weights associated with \eqn{n} observations.}
\item{X}{normalized design matrix. Useful when use the \code{convexity} function in this package.}
\item{betalasso}{The fitted matrix of coefficients with LASSO penalty. The number of rows is equal to the number of coefficients in the model, and the number of columns is the number of \eqn{\lambda}.}
\item{betamcp}{The fitted matrix of coefficients with MCP penalty if \code{penalty = "MCP"}. The number of rows is equal to the number of coefficients in the model, and the number of columns is \code{n_lambda} times \code{n_kappa}. The \eqn{u}th column is the solution when \eqn{\lambda} equals the \eqn{u mod n_lambda}th element of vector \code{lambda} and \eqn{\kappa} equals the \code{floor}(u / n_lambda + 1)th element of vector \code{kappa}. See \code{lambda} and \code{kappa} below for more information.}
\item{betascad}{The fitted matrix of coefficients with SCAD penalty if \code{penalty = "SCAD"}. The number of rows is equal to the number of coefficients in the model, and the number of columns is \code{n_lambda} times \code{n_kappa}. The \eqn{u}th column is the solution when \eqn{\lambda} equals the \eqn{u mod n_lambda}th element of vector \code{lambda} and \eqn{\kappa} equals the \code{floor}(u / n_lambda + 1)th element of vector \code{kappa}. See \code{lambda} and \code{kappa} below for more information.}
\item{betaadap}{The fitted matrix of coefficients with adaptive LASSO penalty if \code{penalty = "adaptive"}. The number of rows is equal to the number of coefficients in the model, and the number of columns is \code{n_lambda} times \code{n_kappa}. The \eqn{u}th column is the solution when \eqn{\lambda} equals the \eqn{u mod n_lambda}th element of vector \code{lambda} and \eqn{\kappa} equals the \code{floor}(u / n_lambda + 1)th element of vector \code{kappa}. See \code{lambda} and \code{kappa} below for more information.}
\item{iter}{number of iterations in each MCP/SCAD solution.}
\item{lambda}{vector of length \code{n_lambda}. The sequence of parameter \eqn{\lambda} in the grid.}
\item{kappa}{vector of length \code{n_kappa}. The sequence of parameter \eqn{\kappa} in the grid.}
\item{n_lambda}{number of lambda values in the grid.}
\item{n_kappa}{number of kappa values in the grid.}
\item{r}{vector of length \eqn{p}. It contains normlizing constants for each column of design matrix \eqn{X}.}
\item{n}{number of rows of matrix \eqn{X}.}
\item{p}{number of columsn of matrix \eqn{X}.}
\item{name}{The names associated with the selected variables.}}
\author{Hao Chai <\email{hao.chai@yale.edu}>}
\arguments{\item{X}{Design matrix of size \eqn{n} by \eqn{p}, where \eqn{n} is the sample size and \eqn{p} is the number of variables.}
\item{y}{response vector of length \eqn{n}.}
\item{delta}{indicator vector of length \eqn{n}. Useful when \code{AFT} is \code{TRUE}. If \eqn{\delta_i} equals 1, \eqn{y_i} is log failure time. If \eqn{\delta_i} equals 0, \eqn{y_i} is log censored time.}
\item{n_lambda}{number of lambda values in the grid. Default value is
  100. Either \code{delta} or \code{weight} has to be specified by
  user.}
\item{n_kappa}{number of kappa values in the grid. Default value is 1.}
\item{kappa0}{Maximum value of \eqn{\kappa}. If \code{penalty = "MCP" or "SCAD"}, it defaults 1/3, and the range of \code{kappa} is (0, kappa0). If \code{penalty = "adaptive"}, it defaults 2, and range of \code{kappa} is (1, kappa0). In the adaptive LASSO case, \code{kappa} is the power to which the initial estimator is raised.}
\item{lambda0}{Maximum value of \eqn{\lambda}. Default value is -1. If \code{lambda0} = -1, the maximum value of \eqn{\lambda} will be the minimum \eqn{\lambda} such that all estimates of \eqn{\beta} are zero.}
\item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than eps. Default is 1e-6.}
\item{max.iter}{Maximum number of iteration. Default is 100.}
\item{method}{Character. Either "exclude" or "include". If \code{method}
  equals "exclude", after the LASSO step, the function will
  automatically exclude all the zero estimates from LASSO and pass only
  the variables associated with non-zero estimates to MCP step. If
  \code{method} equals "include", then LASSO estimates will serve as a
  warm start value and all variables will be included in the MCP
  step. Default is "include".}
\item{name.list}{vector to specify the names of rows of matrix \code{X}. Default is \code{NA}.}
\item{weight}{Only useful when \code{AFT} is \code{FALSE}. Specifies the
  user defined weights associated with observations. Default is
  \code{NA}. Either \code{delta} or \code{weight} has to be specified by
  user.}
\item{AFT}{logical. If \code{TRUE}, the function uses the accelerated
  failure time model to fit the response variable and \code{weight} will
  be calculated as Kaplan-Meier weights. If \code{FALSE}, \code{weight}
  will be the user specified weights for each observation. Default is
  \code{TRUE}.}
\item{penalty}{The penalty to be added to the objective function. Values could be "MCP" (the default), "SCAD" or "adaptive" for MCP, SCAD and adaptive LASSO penalties, respectively. LASSO result is automatically computed.}
}
\examples{X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X \%*\% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta,
n_kappa = 40, kappa0 = 0.99)
}
