paths.plot = function(result, penalty = 'MCP', n_kappa = 1)
{
  col <- hsv(seq(0, 1, len = (result$p + 1)))[1:result$p]
  if ((penalty == "MCP") || (penalty == "mcp"))
  {
    if (n_kappa > result$n_kappa) stop("Regression is not computed for kappa")
    l = result$n_lambda
    init = (n_kappa - 1) * l
    ylim = max(abs(result$betamcp[, (init + 1) : (init + l)]))
    plot(result$betamcp[1, (init + 1) : (init + l)] ~ result$lambda, type = 'l',
         col = col[1], xlim = c(max(result$lambda), 0), 
         ylim = c(-ylim * 1.2, ylim * 1.2), lwd = 1.5,
         main = paste('MCP Path kappa = ', round(result$kappa[n_kappa], 2)), 
         ylab = expression(beta), xlab = expression(lambda))
    for (i in 2 : nrow(result$betamcp))
    lines(result$betamcp[i, (init + 1) : (init + l)] ~ result$lambda, 
          type = 'l', col = col[i], lwd = 1.5)
    grid(10, 10)
  }
  else if ((penalty == "SCAD") || (penalty == "scad"))
  {
    if (n_kappa > result$n_kappa) stop("Regression is not computed for kappa")
    l = result$n_lambda
    init = (n_kappa - 1) * l
    ylim = max(abs(result$betascad[, (init + 1) : (init + l)]))
    plot(result$betascad[1, (init + 1) : (init + l)] ~ result$lambda, type = 'l',
         col = col[1], xlim = c(max(result$lambda), 0), 
         ylim = c(-ylim * 1.2, ylim * 1.2), lwd = 1.5, 
         main = paste('SCAD Path kappa = ', round(result$kappa[n_kappa], 2)), 
         ylab = 'beta', xlab = 'lambda')
    for (i in 2 : nrow(result$betascad))
      lines(result$betascad[i, (init + 1) : (init + l)] ~ result$lambda, 
            type = 'l', col = col[i], lwd = 1.5)
    grid(10, 10)
  }
  else if ((penalty == "LASSO") || (penalty == "lasso"))
  {
    plot(result$betalasso[1, ] ~ result$lambda, type = 'l', col = col[1], 
         xlim = c(max(result$lambda), 0), ylim = c(-max(abs(result$betalasso)) * 
           1.2, max(abs(result$betalasso)) * 1.2), main = 'LASSO Path', 
         ylab = 'beta', xlab = 'lambda', lwd = 1.5)
    for (i in 2:nrow(result$betalasso))
    lines(result$betalasso[i, ] ~ result$lambda, type = 'l', col = col[i],
          lwd = 1.5)
    grid(10, 10)
  }
  else if ((penalty == "adaptive") || (penalty == "ADAPTIVE"))
  {
    if (n_kappa > result$n_kappa) stop("Regression is not computed for kappa")
    l = result$n_lambda
    init = (n_kappa - 1) * l
    ylim = max(abs(result$betaadap[, (init + 1) : (init + l)]))
    plot(result$betaadap[1, (init + 1) : (init + l)] ~ result$lambda, type = 'l',
         col = col[1], xlim = c(max(result$lambda), 0), 
         ylim = c(-ylim * 1.2, ylim * 1.2), lwd = 1.5, 
         main = paste('adaptive-LASSO Path power = ', round(result$kappa[n_kappa], 2)), 
         ylab = 'beta', xlab = 'lambda')
    for (i in 2 : nrow(result$betaadap))
      lines(result$betaadap[i, (init + 1) : (init + l)] ~ result$lambda, 
            type = 'l', col = col[i], lwd = 1.5)
    grid(10, 10)
  }
}
