weighted.sparsenet = function(X, y, delta, n_lambda = 100, n_kappa = 1, 
                              kappa0 = ifelse(penalty == "MCP" || penalty == "SCAD", 
                                              1 / 3, 2), lambda0 = -1, 
                              eps = 1e-6, 
                              max.iter = 1000, method = "include", 
                              name.list = NA, weight = NA, AFT = TRUE, 
                              penalty = "MCP") 
{
  if (AFT) nalist = (is.na(y) || is.na(delta)) else nalist = is.na(y)
  if (sum(nalist) != 0) {
    warning("Missing values in the response have been removed!")
    y = y[!nalist]
    X = X[!nalist, ]
    if (AFT) delta = delta[!nalist] else weight = weight[!nalist]    
  }
  if (sum(is.na(X)) != 0) {
    warning("Missing values in the covariates have been removed!")
    X = X[, apply(X, 2, function(x) !any(is.na(x))), drop = F]
    if (!all(is.na(name.list))) 
      name.list = name.list[apply(X, 2, function(x) !any(is.na(x)))]
  }

  n = nrow(X)
  if (n == 1) stop("n = 1")
  if (n != length(y)) stop("Unmatched dimensionality")
  if (AFT) if (n != length(delta)) stop("Unmatched dimensionality")
  if (kappa0 >= 1 && (penalty == "SCAD" || penalty == "SCAD")) 
    stop("kappa must be strictly less than 1 for MCP or SCAD penalty.")
  if (kappa0 < 1 && penalty == "adaptive") 
    stop("kappa must be greater than or equal to 1 for adaptive LASSO penalty.")
  if (is.vector(X) == T) p = 1 else p = ncol(X)

  orders = sort(y, index.return = T)$ix
  if (AFT) delta = delta[orders]
  y = y[orders]
  X = X[orders, ]
  if (AFT) {
    tmp = ((n - 1) / n) ^ delta[1]
    for (i in 2:(n - 1)) 
      tmp[i] = tmp[i - 1] * ((n - i) / (n - i + 1)) ^ delta[i]
    omega = delta[1] / n
    for (i in 2:n) 
      omega[i] = delta[i] / (n - i + 1) * tmp[i - 1]
    rm(tmp)
  } else {omega = weight[orders]}
  v = sqrt(omega)
  ######## normalization #######

  meany = as.numeric(crossprod(omega, y) / sum(omega))
  meanX = as.numeric(crossprod(omega, X) / sum(omega))
  y = y - meany
  y = v * y
  X = sweep(X, 2, meanX)
  X = sweep(X, 1, v, "*")
  r = sqrt(colSums(X * X))
  X = t(t(X) / r)
  
  if (lambda0 == -1) lambda0 = max(abs(as.vector(y %*% X)))
  beta = matrix(rep(0, n_lambda * p), nrow = n_lambda)
  if (p >= n) 
  {
    lambda = exp(seq(log(lambda0), log(lambda0*0.05), len = n_lambda))
  }  else 
  {lambda = exp(seq(log(lambda0), log(lambda0 * 0.001), len = n_lambda))}
  iter = rep(0, n_lambda * n_kappa)
  if (penalty == "MCP" || penalty == "SCAD") 
    kappa = seq(0, kappa0, len = n_kappa + 1)[2:(n_kappa + 1)]
  if (penalty == "adaptive") 
    kappa = seq(1, kappa0, len = n_kappa + 1)[2:(n_kappa + 1)]
  result = .Call("pathSearch", X, y, method, lambda, kappa, n_lambda, n_kappa, 
                 eps, max.iter, penalty)
  result[[2]] = result[[2]] / r
  result[[1]] = result[[1]] / r
  estimation = list()
  estimation$omega = omega
  estimation$X = X
  estimation$betalasso = result[[1]]
  if (penalty == "MCP") estimation$betamcp = result[[2]] else
    if (penalty == "SCAD") estimation$betascad = result[[2]] else
      if (penalty == "adaptive") estimation$betaadap = result[[2]]
  estimation$iter = result[[3]]
  estimation$lambda = result[[4]]
  estimation$kappa = result[[5]]
  estimation$n_lambda = n_lambda
  estimation$n_kappa = n_kappa
  estimation$r = r
  estimation$n = n
  estimation$p = p
  if (!all(is.na(name.list))) estimation$name = name.list
  return(estimation)
}