bic.sparsenet = function(result, X, y)
{
  n = result$n
  p = result$p
  losslasso = log(result$omega %*% 
    (sweep(crossprod(t(X), result$betalasso), 1, y, "-")) ^ 2) + 
    apply((result$betalasso != 0), 2, sum) * log(n) / n
  indexlasso = which.min(losslasso)
  lossmcp = log(result$omega %*% 
    (sweep(crossprod(t(X), result$betamcp), 1, y, "-")) ^ 2) + 
    apply((result$betamcp != 0), 2, sum) * log(n) / n
  indexmcp = which.min(lossmcp)
  best.kappa = floor((indexmcp - 1) / result$n_lambda) + 1
  best.lambda.mcp = result$lambda[indexmcp - (best.kappa - 1) * result$n_lambda]
  best.kappa = result$kappa[best.kappa]
  return(list(betamcp = result$betamcp[, indexmcp], 
              betalasso = result$betalasso[, indexlasso], 
              lasso.lambda = result$lambda[indexlasso], 
              mcp.lambda = best.lambda.mcp, mcp.kappa = best.kappa, 
              name.lasso = result$name[abs(result$betalasso[, indexlasso]) >= 
                1e-06], name.mcp = result$name[abs(result$betamcp[, indexmcp]) 
                                               >= 1e-06]))
}