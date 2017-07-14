cv.sparsenet = function(result, X, y,  delta, n_lambda = 100, n_kappa = 1, 
                        kappa0 = 1 / 3, lambda0 = -1, eps = 1e-6, 
                        max.iter = 100, method = "include", power = 2, 
                        fold = 5, weight = NA, AFT = TRUE, penalty = "MCP", 
                        stratified = TRUE)
{
  if (AFT) nalist = (is.na(y) || is.na(delta)) else nalist = is.na(y)
  if (sum(nalist) != 0) {
    warning("Missing values in the response have been removed!")
    y = y[!nalist]
    X = X[!nalist,]
    if (AFT) delta = delta[!nalist] else weight = weight[!nalist]
  }
  if (sum(is.na(X)) != 0) {
    warning("Missing values in the covariates have been removed!")
    X = X[, apply(X, 2, function(x) !any(is.na(x))), drop = F]
  }
  sum.resid.mcp = rep(0, n_lambda * n_kappa)
  sum.resid.lasso = rep(0, n_lambda)
  N = nrow(X)
  if (((sum(delta) < 2 * fold) | (N - sum(delta) < 2 * fold)) & stratified)
  {
    stratified = F
    warning("Stratified sampling is disabled because one group is too small.")
  }
  if (stratified)
  {
    size = floor(N / fold)
    fail.size = round(size * mean(delta))
    cens.size = size - fail.size
    if (cens.size == 0) {cens.size = 1; fail.size = size - cens.size}
    if (fail.size == 0) {fail.size = 1; cens.size = size - fail.size}
    fail.index = sample(seq(1, N)[delta == 1])
    cens.index = sample(seq(1, N)[delta == 0])
  } else
  {
    random.index = sample(seq(1, N), replace = F)
    size = floor(N / fold)
  }
  for (i in 1:(fold - 1)){
    if (stratified)
    {
      valid.index = c(fail.index[((i - 1) * fail.size + 1):(i * fail.size)], 
                      cens.index[((i - 1) * cens.size + 1):(i * cens.size)])
      train.index = c(fail.index[-(((i - 1) * fail.size + 1):(i * fail.size))], 
                      cens.index[-(((i - 1) * cens.size + 1):(i * cens.size))])
    } else
    {
      valid.index = random.index[((i - 1) * size + 1):(i * size)]
      train.index = random.index[-(((i - 1) * size + 1):(i * size))]
    }
    yy = y[valid.index]
    XX = X[valid.index, ]
    ddelta = delta[valid.index]
    orders = sort(yy, index.return = T)$ix
    ddelta = ddelta[orders]
    yy = yy[orders]
    XX = XX[orders, ]
    nn = size
    if (AFT) {
      tmp = ((nn - 1) / nn) ^ ddelta[1]
      for (j in 2:(nn - 1)) 
        tmp[j] = tmp[j - 1] * ((nn - j) / (nn - j + 1)) ^ ddelta[j]
      omega = delta[1] / nn
      for (j in 2:nn) 
        omega[j] = ddelta[i] / (nn - j + 1) * tmp[j - 1]
      rm(tmp)
    } else
    {
      wweight = weight[valid.index]
      omega = wweight[orders]
    }
    if (AFT) {fit = weighted.sparsenet(X = X[train.index, ], y = y[train.index], 
                                       delta = delta[train.index], n_lambda = n_lambda, 
                                       n_kappa = n_kappa, kappa0 = kappa0, lambda0 = -1, eps = eps, 
                                       max.iter = max.iter, method = method, AFT = AFT, 
                                       penalty = penalty)} else
                                       {fit = weighted.sparsenet(X = X[train.index, ], y = y[train.index], 
                                                                 n_lambda = n_lambda, n_kappa = n_kappa, kappa0 = kappa0, 
                                                                 lambda0 = -1, eps = eps, max.iter = max.iter, method = method, 
                                                                 weight = weight[train.index], AFT = F, penalty = penalty
                                                                 )}

    meanX = as.numeric(crossprod(fit$omega, X[train.index, ]) / sum(fit$omega))
    meany = as.numeric(crossprod(fit$omega, y[train.index])) / sum(fit$omega)
    if (penalty == "MCP") 
    {
      betamcp0 = as.vector(meany - crossprod(meanX, fit$betamcp))
      betamcp0 = matrix(rep(betamcp0, length(valid.index)), 
                        ncol = length(betamcp0), byrow = T)
      predicted.mcp = crossprod(t(XX), fit$betamcp) + betamcp0
      sum.resid.mcp = sum.resid.mcp + crossprod(omega, abs(predicted.mcp - yy) 
                                                ^ power)
    }
    if (penalty == "SCAD")
    {
      betascad0 = as.vector(meany - crossprod(meanX, fit$betascad))
      betascad0 = matrix(rep(betascad0, length(valid.index)), 
                         ncol = length(betascad0), byrow = T)
      predicted.scad = crossprod(t(XX), fit$betascad) + betascad0
      sum.resid.scad = sum.resid.scad + 
        crossprod(omega, abs(predicted.scad - yy) ^ power)
    }
    if (penalty == "adaptive")
    {
      betaadap0 = as.vector(meany - crossprod(meanX, fit$betaadap))
      betaadap0 = matrix(rep(betaadap0, length(valid.index)), 
                         ncol = length(betaadap0), byrow = T)
      predicted.adap = crossprod(t(XX), fit$betaadap) + betaadap0
      sum.resid.adap = sum.resid.adap + 
        crossprod(omega, abs(predicted.adap - yy) ^ power)
    }
    betalasso0 = as.vector(meany - crossprod(meanX, fit$betalasso))
    betalasso0 = matrix(rep(betalasso0, length(valid.index)), ncol = length(betalasso0), 
                        byrow = T)
    predicted.lasso = crossprod(t(XX), fit$betalasso) + betalasso0
    sum.resid.lasso = sum.resid.lasso + crossprod(omega, abs(predicted.lasso - yy) ^ 
                                                    power)
    rm(omega, orders, yy, XX, ddelta, nn)
  }
  if (stratified)
  {
    valid.index = c(fail.index[((fold - 1) * fail.size + 1):length(fail.index)], 
                    cens.index[((fold - 1) * cens.size + 1):length(cens.index)])
    train.index = c(fail.index[-(((fold - 1) * fail.size + 1):length(fail.index))], 
                    cens.index[-(((fold - 1) * cens.size + 1):length(cens.index))])
  } else
  {
    valid.index = random.index[((fold-1) * size + 1):N]
    train.index = random.index[-(((fold-1) * size + 1):N)]
  }
  yy = y[valid.index]
  XX = X[valid.index, ]
  ddelta = delta[valid.index]
  orders = sort(yy, index.return = T)$ix
  ddelta = ddelta[orders]
  yy = yy[orders]
  XX = XX[orders, ]
  nn = N - (fold - 1) * size
  if (AFT) {
    tmp = ((nn - 1) / nn) ^ ddelta[1]
    for (j in 2:(nn - 1)) 
      tmp[j] = tmp[j - 1] * ((nn - j) / (nn - j + 1)) ^ ddelta[j]
    omega = delta[1] / nn
    for (j in 2:nn) 
      omega[j] = ddelta[j] / (nn - j + 1) * tmp[j - 1]
    rm(tmp)
  } else
  {
    wweight = weight[valid.index]
    omega = wweight[orders]
  }
  if (AFT) {fit = weighted.sparsenet(X = X[train.index, ], y = y[train.index], 
                                     delta = delta[train.index], n_lambda = n_lambda, 
                                     n_kappa = n_kappa, kappa0 = kappa0, lambda0 = -1, eps = eps, 
                                     max.iter = max.iter, method = method, AFT = AFT, 
                                     penalty = penalty)} else
                                     {fit = weighted.sparsenet(X = X[train.index, ], y = y[train.index], 
                                                               n_lambda = n_lambda, n_kappa = n_kappa, kappa0 = kappa0, 
                                                               lambda0 = -1, eps = eps, max.iter = max.iter, method = method, 
                                                               weight = weight[train.index], AFT = F, penalty = penalty)}

  meanX = as.numeric(crossprod(fit$omega, X[train.index, ]) / sum(fit$omega))
  meany = as.numeric(crossprod(fit$omega, y[train.index])) / sum(fit$omega)
  if (penalty == "MCP") 
  {
    betamcp0 = as.vector(meany-crossprod(meanX, fit$betamcp))
    betamcp0 = matrix(rep(betamcp0, length(valid.index)), 
                      ncol = length(betamcp0), byrow = T)
    predicted.mcp = crossprod(t(XX), fit$betamcp) + betamcp0
    sum.resid.mcp = sum.resid.mcp + crossprod(omega, abs(predicted.mcp - yy) ^ power)
  }
  if (penalty == "SCAD")
  {
    betascad0 = as.vector(meany - crossprod(meanX, fit$betascad))
    betascad0 = matrix(rep(betascad0, length(valid.index)), 
                       ncol = length(betascad0), 
                       byrow = T)
    predicted.scad = crossprod(t(XX), fit$betascad) + betascad0
    sum.resid.scad = sum.resid.scad + 
      crossprod(omega, abs(predicted.scad - yy) ^ power)
  }
  if (penalty == "adaptive")
  {
    betaadap0 = as.vector(meany - crossprod(meanX, fit$betaadap))
    betaadap0 = matrix(rep(betaadap0, length(valid.index)), 
                       ncol = length(betaadap0), 
                       byrow = T)
    predicted.adap = crossprod(t(XX), fit$betaadap) + betaadap0
    sum.resid.adap = sum.resid.adap + 
      crossprod(omega, abs(predicted.adap - yy) ^ power)
  }
  betalasso0 = as.vector(meany - crossprod(meanX, fit$betalasso))
  betalasso0 = matrix(rep(betalasso0, length(valid.index)), 
                      ncol = length(betalasso0), 
                      byrow = T)
  predicted.lasso = crossprod(t(XX), fit$betalasso) + betalasso0
  sum.resid.lasso = sum.resid.lasso + crossprod(omega, abs(predicted.lasso - yy) ^ 
                                                  power)
  rm(omega, orders, yy, XX, ddelta, nn)
  if (penalty == "MCP") 
  {
    indexmcp = which.min(sum.resid.mcp)
    best.kappa = floor((indexmcp - 1) / result$n_lambda) + 1
    best.lambda = result$lambda[indexmcp - (best.kappa - 1) * result$n_lambda]
  }
  if (penalty == "SCAD") 
  {
    indexscad = which.min(sum.resid.scad)
    best.kappa = floor((indexscad -1) / result$n_lambda) + 1
    best.lambda = result$lambda[indexscad - (best.kappa - 1) * result$n_lambda] 
  }
  if (penalty == "adaptive") 
  {
    indexadap = which.min(sum.resid.adap)
    best.kappa = floor((indexadap -1) / result$n_lambda) + 1
    best.lambda = result$lambda[indexadap - (best.kappa - 1) * result$n_lambda] 
  }
  best.kappa = result$kappa[best.kappa]
  indexlasso = which.min(sum.resid.lasso)
  lasso.lambda = result$lambda[indexlasso]
  if (penalty == "MCP")
    return.list = list(betamcp = result$betamcp[, indexmcp], 
                       betalasso = result$betalasso[, indexlasso], 
                       mcp.lambda = best.lambda, mcp.kappa = best.kappa, 
                       lasso.lambda = lasso.lambda, 
                       name.lasso = 
                         result$name[abs(result$betalasso[, indexlasso]) >= 1e-06], 
                       name.mcp = result$name[abs(result$betamcp[, indexmcp]) >= 1e-06])
  if (penalty == "SCAD")
    return.list = list(betascad = result$betascad[, indexscad], 
                       betalasso = result$betalasso[, indexlasso], 
                       scad.lambda = best.lambda, scad.kappa = best.kappa, 
                       lasso.lambda = lasso.lambda, 
                       name.lasso = 
                         result$name[abs(result$betalasso[, indexlasso]) >= 1e-06], 
                       name.scad = 
                         result$name[abs(result$betascad[, indexscad]) >= 1e-06])
  if (penalty == "adaptive")
    return.list = list(betaadap = result$betaadap[, indexadap], 
                       betalasso = result$betalasso[, indexlasso], 
                       adap.lambda = best.lambda, adap.kappa = best.kappa, 
                       lasso.lambda = lasso.lambda, 
                       name.lasso = 
                         result$name[abs(result$betalasso[, indexlasso]) >= 1e-06], 
                       name.adap = 
                         result$name[abs(result$betaadap[, indexadap]) >= 1e-06])
  return(return.list)
}