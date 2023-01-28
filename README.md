# MZOIBTS

Beta(mu,phi): This function includes four sub-functions: dbeta0, pbeta0, qbeta0, and rbeta0 that correspond to the density, cumulative distribution, quantile and random number generation functions for the beta distribution, respectively. These functions take in x, mu, and phi as arguments.

ZOIB(p1,p2,mu,phi): This function includes four sub-functions: dZOIB, pZOIB, qZOIB, and rZOIB that correspond to the density, cumulative distribution, quantile and random number generation functions for the ZOIB distribution, respectively. These functions take in x, p1, p2, mu, and phi as arguments.

rZOIBTS: This function generates a Zero-One-inflated Beta Time Series (ZOIBTS) with n observations. It takes in n, copula, alpha, p1, p2, mu, and phi as arguments.

rMZOIBTS: This function generates a Marginalized Zero-One-inflated Beta Time Series (MZOIBTS) with n observations. It takes in n, copula, alpha, X1, beta1, X2, beta2, X3, beta3, X4, and beta4 as arguments.

est.ILL: This function outputs the independence-loglikelihood estimator of a Marginalized Zero-One-inflated Beta Time Series (MZOIBTS) model. It takes in X1, X2, X3, X4, and y as arguments.

qest.cop: This function outputs the quasi-likelihood estimator of the copula parameter. It takes in theta (outputs from est.ILL), X1, X2, X3, X4, y, and copula as arguments. Available copula families are Gaussian, Clayton, Gumbel, Frank, and AMH.

twostep.MZOIBTS: This function fits a Marginalized Zero-One-inflated Beta Time Series (MZOIBTS) model. It takes in R (default 100), X1, X2, X3, X4, y, and copula as arguments. Available copula families are Gaussian, Clayton, Gumbel, Frank, and AMH.

selection: This selects the optimal change points for a Marginalized Zero-One-inflated Beta Interupted Time Series (MZOIBITS) model. It takes in X1, X2, y, and candidates as arguments.
