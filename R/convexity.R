convexity = function(result, eps = sqrt(.Machine$double.eps), 
                     eigen.min = 1e-3, concave.rm = TRUE)
{
  return.list = list(lambda = result$lambda, kappa = result$kappa, 
                     index = matrix())
  index = matrix(rep(0, 4 * result$n_lambda * result$n_kappa), ncol = 4)
  colnames(index)= c("lambda.index", "kappa.index", "beta.index", "min.eigen")
  l = result$n_lambda
  subset = seq(1:l)
  xProd = crossprod(result$X, result$X)
  r = result$r
  count = 1
  for (i in 1:result$n_kappa)
  {
    tmp = result$betamcp[, seq((i - 1) * result$n_lambda + 1, 
                               length.out = result$n_lambda)]
    tmp = tmp[, subset]
    indMatrix = (abs(tmp) > eps)
    J = length(subset)
    kappa_ = result$kappa[i]
    if (concave.rm) subset.ind = NULL
    eigen.value = rep(NA, J)
    for (j in 1:J)
    {
      lambda_ = result$lambda[subset[j]]
      beta = tmp[indMatrix[, j], j]
      secondD = rep(0, sum(indMatrix[, j]))
      secondD[abs(beta) < (lambda_ * kappa_ / r[indMatrix[, j]])] = kappa_ * 
        if (length(secondD) <= 1) 
        {
          minEigen = as.numeric((xProd[indMatrix[, j], indMatrix[, j]] - 
            secondD) * r[indMatrix[, j]] ^ 2)
        } else
        {
          minEigen = min(eigen(t(xProd[indMatrix[, j], indMatrix[, j]] - 
            diag(secondD)) * r[indMatrix[, j]] ^ 2, only.values = T, 
                               symmetric = T)$values)
        }
      if (length(minEigen) > 0)
      {
        if ((minEigen < eigen.min) && concave.rm) subset.ind = c(subset.ind, j)
        eigen.value[j] = minEigen
      }
    }
    rm(tmp)
    rm(indMatrix)
    if (concave.rm){
      if (length(subset.ind) > 0) 
      {
        subset = subset[-subset.ind]
        eigen.value = eigen.value[-subset.ind]
      }
      if (length(subset) == 0) break
    }
    for (k in 1 : length(subset))
    {
      index[count, 1] = subset[k]
      index[count, 2] = i
      index[count, 3] = result$n_lambda * (i - 1) + subset[k]
      index[count, 4] = eigen.value[k]
      count = count + 1
    }
  }
  count = count - 1
  index = index[1:count, ]
  return.list$index = index
  return.list
}