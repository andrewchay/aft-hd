library(aftHD)
X = matrix(rnorm(8000), nrow = 20)
beta0 = c(rep(10, 5), rep(0, 395))
y = rnorm(20) + X %*% beta0
delta = rep(1, 20)
result = weighted.sparsenet(X, y, delta, n_kappa = 40, kappa0 = 0.99)
bic.result = bic.sparsenet(result, X, y)
cv.result = cv.sparsenet(result, X, y, delta, n_kappa = 40, kappa0 = 0.99)
a = convexity(result, concave.rm = F)
# library(spatstat)
# u = a$index[, c(1, 2)]
# uu = matrix(rep(0, result$n_lambda * result$n_kappa), nrow = 100)
# uu[u] = a$index[, 4]
# im.uu = im(uu)
# plot(im.uu, ylim = 0:100, xlab = expression(kappa), 
#      ylab = expression(lambda), col = heat.colors(10), 
#      main = "Minimal eigen value over the tuning parameter grid")
# axis(1, c(1, result$n_kappa))
# axis(2, c(1, result$n_lambda), pos = -10)
# paths.plot(result, n_kappa = floor(result$n_kappa / 2))
# paths.plot(result, "lasso")